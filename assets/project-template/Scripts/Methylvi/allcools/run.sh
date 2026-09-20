#!/usr/bin/env bash
# Unified .cov -> ALLCools -> MethylVI entry point.
set -euo pipefail
# The run DAG is the production path and is identified by SCMO_MANAGED_RUN_ID,
# which task_adapter.py injects. A standalone run is allowed only when it is
# explicitly acknowledged, so a non-DAG run is never reached by accident. That
# mode writes no task_outputs.json and is not a run's completion evidence; it
# exists for exploring a stage, smoke testing a change, and recovering one stage.
if [[ -z "${SCMO_MANAGED_RUN_ID:-}" ]]; then
  if [[ "${SCMO_STANDALONE_ACK:-0}" != 1 ]]; then
    echo "Production runs go through the run DAG: plan_workflow.py, then submit_workflow.py." >&2
    echo "For exploration, smoke testing, or recovery, set SCMO_STANDALONE_ACK=1." >&2
    exit 2
  fi
  echo "WARNING: standalone run without a managed run ID; no completion evidence is recorded." >&2
fi
# The acknowledgement above says a run may bypass the DAG; it says nothing about
# which host may execute it. A submit host has a controller to submit to, so the
# work belongs in a job there. This is the shell counterpart of the guard in
# run_task.py, and it exists because this script is reachable by hand: without it
# a direct call runs the whole route on the login node. An allocation sets
# SLURM_JOB_ID, so the packaged sbatch wrappers pass through untouched.
authorised=0
if [[ -n "${SCMO_MANAGED_RUN_ID:-}" && "${SCMO_LOGIN_EXECUTION_ACK:-}" == "${SCMO_MANAGED_RUN_ID}" ]]; then
  authorised=1
fi
if [[ -z "${SLURM_JOB_ID:-}" && "$authorised" == 0 ]] && ping=$(scontrol ping 2>/dev/null) && [[ "$ping" == *"is UP"* ]]; then
  echo "Refusing to run: $(hostname) is a Slurm submit host (login node), so this route would" >&2
  echo "run on the login node. Submit the run DAG instead: submit_workflow.py, or work inside an" >&2
  echo "allocation: salloc, or srun --pty bash." >&2
  exit 2
fi

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# The run DAG and the packaged sbatch wrappers export SCMO_PROJECT_ROOT; a call
# made by hand does not, so locate the project from this script's own path. This
# mirrors the submit_*.sh entry points and keeps the direct call described in
# references/execution.md working without extra setup.
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$here/../../.." && pwd)}
source "$here/00_methylvi_config.sh"
shared_scripts="$here/../shared"
stage=${1:-all}

# ALLCools validates helper programs (tabix/bgzip/bedtools) through PATH.
export PATH="$SCMO_ALLCOOLS_ENV/bin:$PATH"
numeric_threads="${SCMO_NUMERIC_THREADS:-1}"
export OMP_NUM_THREADS="$numeric_threads"
export MKL_NUM_THREADS="$numeric_threads"
export OPENBLAS_NUM_THREADS="$numeric_threads"
export NUMEXPR_NUM_THREADS="$numeric_threads"
export MPLBACKEND="${MPLBACKEND:-Agg}"

allc_python="$SCMO_ALLCOOLS_ENV/bin/python"
allcools="$SCMO_ALLCOOLS_ENV/bin/allcools"
mvi_python="$SCMO_METHYLVI_ENV/bin/python"

require_file() {
  [[ -s "$1" ]] || { echo "ERROR: required file is missing or empty: $1" >&2; exit 1; }
}

verify() {
  [[ -d "$SCMO_COV_DIR" ]] || { echo "ERROR: coverage directory missing: $SCMO_COV_DIR" >&2; exit 1; }
  require_file "$SCMO_ANNOTATION"
  require_file "$SCMO_CHROM_SIZES"
  require_file "$SCMO_BLACKLIST"
  [[ -x "$allc_python" && -x "$allcools" ]] || { echo "ERROR: invalid ALLCools environment: $SCMO_ALLCOOLS_ENV" >&2; exit 1; }
  [[ -x "$mvi_python" ]] || { echo "ERROR: invalid MethylVI environment: $SCMO_METHYLVI_ENV" >&2; exit 1; }
  command -v bedtools >/dev/null || { echo "ERROR: bedtools is required for blacklist filtering" >&2; exit 1; }
  command -v intersectBed >/dev/null || { echo "ERROR: intersectBed is required for blacklist filtering" >&2; exit 1; }
  "$allc_python" -c 'import ALLCools,anndata,pandas,scanpy; print("ALLCools environment OK", ALLCools.__version__)'
  "$mvi_python" -c 'import anndata,mudata,scanpy,scvi,torch; from scvi.external import METHYLVI; print("MethylVI environment OK", scvi.__version__)'
  "$allc_python" "$here/01_prepare_allcools.py" --verify-only
}

prepare_allc() {
  "$allc_python" "$here/01_prepare_allcools.py"
}

generate_mcds() {
  require_file "$SCMO_ALLC_TABLE"
  mkdir -p "$SCMO_ALLCOOLS_ROOT"
  table_hash=$(sha256sum "$SCMO_ALLC_TABLE" | awk '{print $1}')
  config="table_sha256=$table_hash bin_size=$SCMO_BIN_SIZE context=$SCMO_MC_CONTEXT hypo_cutoff=$SCMO_HYPO_SCORE_CUTOFF"
  config_file="$SCMO_ALLCOOLS_ROOT/mcds.config.txt"
  complete="$SCMO_ALLCOOLS_ROOT/mcds.COMPLETE"
  if [[ -e "$complete" ]]; then
    require_file "$config_file"
    [[ $(<"$config_file") == "$config" ]] || {
      echo "ERROR: existing MCDS was built with different inputs/parameters; use a new SCMO_MVI_ROOT" >&2
      exit 1
    }
    echo "Existing compatible MCDS detected: $SCMO_MCDS"
    return
  fi
  if [[ -e "$SCMO_MCDS" ]]; then
    echo "ERROR: incomplete MCDS path exists without completion marker: $SCMO_MCDS" >&2
    exit 1
  fi
  "$allcools" generate-dataset \
    --allc_table "$SCMO_ALLC_TABLE" \
    --output_path "$SCMO_MCDS" \
    --chrom_size_path "$SCMO_CHROM_SIZES" \
    --obs_dim cell \
    --cpu "$SCMO_THREADS" \
    --chunk_size 10 \
    --regions chrom5k "$SCMO_BIN_SIZE" \
    --quantifiers chrom5k hypo-score "$SCMO_MC_CONTEXT" "cutoff=$SCMO_HYPO_SCORE_CUTOFF"
  printf '%s\n' "$config" > "$config_file"
  touch "$complete"
}

cluster_allcools() {
  require_file "$SCMO_ALLCOOLS_ROOT/mcds.config.txt"
  annotation_hash=$(sha256sum "$SCMO_ANNOTATION" | awk '{print $1}')
  mcds_config_hash=$(sha256sum "$SCMO_ALLCOOLS_ROOT/mcds.config.txt" | awk '{print $1}')
  blacklist_hash=$(md5sum "$SCMO_BLACKLIST" | awk '{print $1}')
  cluster_config="mcds_config_sha256=$mcds_config_hash annotation_sha256=$annotation_hash blacklist_md5=$blacklist_hash blacklist_fraction=$SCMO_BLACKLIST_FRACTION bin_cutoff=$SCMO_BINARIZE_CUTOFF target_features=$SCMO_TARGET_FEATURES lsi_components=$SCMO_LSI_COMPONENTS p_cutoff=$SCMO_LSI_P_CUTOFF neighbors=$SCMO_ALLCOOLS_NEIGHBORS leiden=$SCMO_ALLCOOLS_LEIDEN_RESOLUTION repeats=$SCMO_CONSENSUS_LEIDEN_REPEATS consensus_leiden=$SCMO_CONSENSUS_LEIDEN_RESOLUTION seed=$SCMO_SEED"
  cluster_config_file="$SCMO_ALLCOOLS_ROOT/cluster.config.txt"
  if [[ -s "$SCMO_ALLCOOLS_H5AD" ]]; then
    require_file "$cluster_config_file"
    [[ $(<"$cluster_config_file") == "$cluster_config" ]] || {
      echo "ERROR: existing clustered H5AD uses different parameters; use a new SCMO_MVI_ROOT" >&2
      exit 1
    }
    echo "Existing compatible ALLCools H5AD detected: $SCMO_ALLCOOLS_H5AD"
    return
  fi
  "$allc_python" "$here/02_cluster_allcools.py"
  printf '%s\n' "$cluster_config" > "$cluster_config_file"
}

run_allcools() {
  prepare_allc
  generate_mcds
  cluster_allcools
}

build_methylvi() {
  require_file "$SCMO_ALLCOOLS_H5AD"
  "$mvi_python" "$here/03_build_methylvi_from_allcools.py"
}

train_methylvi() {
  require_file "$SCMO_MVI_INPUT"
  "$mvi_python" "$shared_scripts/04_train_methylvi.py"
}

supervised_umap() {
  require_file "$SCMO_MVI_RESULTS/methylvi_embedding.h5ad"
  "$mvi_python" "$shared_scripts/05_plot_supervised_umap.py"
}

plots_before_methylvi() {
  require_file "$SCMO_ALLCOOLS_H5AD"
  "$allc_python" "$here/06_plot_allcools_umap.py"
}

plots_after_methylvi() {
  require_file "$SCMO_MVI_RESULTS/methylvi_embedding.h5ad"
  "$mvi_python" "$here/07_plot_methylvi_umap.py"
}

case "$stage" in
  verify) verify ;;
  prepare) verify; prepare_allc; generate_mcds ;;
  cluster) cluster_allcools ;;
  allcools) verify; run_allcools ;;
  build) build_methylvi ;;
  train) train_methylvi ;;
  plots-before) plots_before_methylvi ;;
  plots-after) plots_after_methylvi ;;
  supervised) supervised_umap ;;
  all) verify; run_allcools; build_methylvi; train_methylvi; supervised_umap ;;
  *) echo "Usage: bash $0 {verify|prepare|cluster|allcools|build|train|plots-before|plots-after|supervised|all}" >&2; exit 2 ;;
esac
