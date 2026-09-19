#!/usr/bin/env bash
# Complete VMR models whose training finished before the plotting-wrapper fix.
# This is intentionally a post-training recovery stage; it never retrains.
# Usage: sbatch [resource options] recover_vmr_outputs.sbatch
# Resources and log paths are injected by the submitter, not declared here.
#SBATCH --job-name=scmo_mvi_vmr_recover

set -Eeuo pipefail
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi
source "$script_dir/00_experiment_config.sh"
mvi_root="$SCMO_MVI_EXPERIMENT_ROOT"
shared_scripts="$script_dir/shared"
python_exe="$SCMO_METHYLVI_ENV/bin/python"
export SCMO_BATCH_KEY="${SCMO_BATCH_KEY:-sample_id}"
export SCMO_CELLTYPE_KEY="${SCMO_CELLTYPE_KEY:-cell_type}"
export SCMO_SEED="${SCMO_SEED:-0}"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-${SCMO_THREADS:-16}}"
export SCMO_SUPERVISED_TARGET_KEY="${SCMO_SUPERVISED_TARGET_KEY:-$SCMO_CELLTYPE_KEY}"
export SCMO_SUPERVISED_TARGET_WEIGHTS="${SCMO_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export SCMO_SUPERVISED_NEIGHBORS="${SCMO_SUPERVISED_NEIGHBORS:-15}"
export SCMO_SUPERVISED_MIN_DIST="${SCMO_SUPERVISED_MIN_DIST:-0.5}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"

for version in $SCMO_VMR_THRESHOLDS; do
  for features in $SCMO_FEATURE_TARGETS; do
    route="$mvi_root/vmr/var_${version}/features_${features}"
    input="$route/input.h5mu"
    results="$route/results"
    [[ -s "$input" && -s "$results/methylvi_embedding.h5ad" && -s "$results/run_summary.json" ]] || {
      echo "Missing completed training output: $route" >&2
      exit 1
    }
    export SCMO_MVI_INPUT="$input" SCMO_MVI_RESULTS="$results" VMR_MVI_RESULTS="$results"
    "$python_exe" "$script_dir/vmr/03_plot_vmr_methylvi_umap.py" \
      --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
    "$python_exe" "$shared_scripts/05_plot_supervised_umap.py"
    "$python_exe" "$shared_scripts/07_plot_methylation_qc.py" \
      --input "$input" --output-dir "$results/qc"
    touch "$route/model.COMPLETE"
  done
done
printf 'Recovered VMR post-training outputs under %s\n' "$mvi_root"
