#!/usr/bin/env bash
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

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# The run DAG and the packaged sbatch wrappers export SCMO_PROJECT_ROOT; a call
# made by hand does not, so locate the project from this script's own path. This
# mirrors the submit_*.sh entry points and keeps the direct call described in
# references/execution.md working without extra setup.
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$here/../../.." && pwd)}
source "$here/00_vmr_methylvi_config.sh"
stage=${1:-all}
python_exe="$VMR_METHYLVI_ENV/bin/python"
shared_scripts="$here/../shared"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export OMP_NUM_THREADS="${VMR_NUMERIC_THREADS:-1}"
export MKL_NUM_THREADS="${VMR_NUMERIC_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${VMR_NUMERIC_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${VMR_NUMERIC_THREADS:-1}"

require_file() {
  [[ -s "$1" ]] || { echo "ERROR: required file missing or empty: $1" >&2; exit 1; }
}

verify() {
  if [[ -n "$VMR_INPUT_MANIFEST" ]]; then
    require_file "$VMR_METHSCAN_RUN_DIR/smooth.COMPLETE"
    require_file "$VMR_FILTERED_CELL_IDS"
  fi
  [[ -n "$VMR_SOURCE_BED" && -s "$VMR_SOURCE_BED" ]] || {
    echo "ERROR: MethSCAn VMR BED is not ready: $VMR_SOURCE_BED" >&2
    exit 1
  }
  if [[ -n "$VMR_INPUT_MANIFEST" ]]; then
    [[ -s "$VMR_INPUT_MANIFEST" ]] || {
      echo "ERROR: MethSCAn selected-ALLC manifest is missing: $VMR_INPUT_MANIFEST" >&2
      exit 1
    }
  else
    [[ -d "$VMR_COV_DIR" && -d "$VMR_EXISTING_ALLC_DIR" ]] || {
      echo "ERROR: legacy coverage or ALLC directory missing" >&2; exit 1;
    }
  fi
  for path in "$VMR_SOURCE_BED" "$VMR_CHROM_SIZES" "$VMR_BLACKLIST" "$VMR_ANNOTATION"; do
    require_file "$path"
  done
  [[ -x "$python_exe" ]] || { echo "ERROR: invalid MethylVI environment: $VMR_METHYLVI_ENV" >&2; exit 1; }
  "$python_exe" -c 'import anndata,mudata,numpy,pandas,scanpy,scvi,torch; from scvi.external import METHYLVI; print("VMR-MethylVI environment OK", scvi.__version__)'
}

prepare() {
  "$python_exe" "$here/01_prepare_vmr_inputs.py"
}

build() {
  require_file "$VMR_FILTERED_BED"
  require_file "$VMR_ALLC_TABLE"
  "$python_exe" "$here/02_build_vmr_methylvi_input.py"
}

train() {
  require_file "$VMR_MVI_INPUT"
  "$python_exe" "$shared_scripts/04_train_methylvi.py"
}

plots() {
  require_file "$VMR_MVI_RESULTS/methylvi_embedding.h5ad"
  "$python_exe" "$here/03_plot_vmr_methylvi_umap.py"
}

supervised() {
  require_file "$VMR_MVI_RESULTS/methylvi_embedding.h5ad"
  "$python_exe" "$shared_scripts/05_plot_supervised_umap.py"
}

case "$stage" in
  verify) verify ;;
  prepare) verify; prepare ;;
  build) build ;;
  train) train ;;
  plots) plots ;;
  supervised) supervised ;;
  all) verify; prepare; build; train; plots; supervised ;;
  *) echo "Usage: bash $0 {verify|prepare|build|train|plots|supervised|all}" >&2; exit 2 ;;
esac
