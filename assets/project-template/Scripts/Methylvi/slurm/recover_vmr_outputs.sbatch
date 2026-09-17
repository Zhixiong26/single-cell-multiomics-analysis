#!/usr/bin/env bash
# Complete VMR models whose training finished before the plotting-wrapper fix.
# This is intentionally a post-training recovery stage; it never retrains.
#SBATCH --ntasks=1

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
root=${SCMO_MVI_EXPERIMENT_ROOT:?set SCMO_MVI_EXPERIMENT_ROOT}
python_exe=${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python
export SCMO_BATCH_KEY=sample_id SCMO_CELLTYPE_KEY=cell_type SCMO_SEED=0
export SCMO_THREADS=45 SCMO_SUPERVISED_TARGET_KEY=cell_type
export SCMO_SUPERVISED_TARGET_WEIGHTS="0.2 0.5 0.7 0.9"
export SCMO_SUPERVISED_NEIGHBORS=15 SCMO_SUPERVISED_MIN_DIST=0.5 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

for version in 0.01 0.02 0.05; do
  for features in 10000 30000; do
    route="$root/vmr/var_${version}/features_${features}"
    input="$route/input.h5mu"
    results="$route/results"
    [[ -s "$input" && -s "$results/methylvi_embedding.h5ad" && -s "$results/run_summary.json" ]] || {
      echo "Missing completed training output: $route" >&2
      exit 1
    }
    export SCMO_MVI_INPUT="$input" SCMO_MVI_RESULTS="$results" VMR_MVI_RESULTS="$results"
    "$python_exe" Scripts/Methylvi/vmr/03_plot_vmr_methylvi_umap.py \
      --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
    "$python_exe" Scripts/Methylvi/shared/05_plot_supervised_umap.py
    "$python_exe" Scripts/Methylvi/shared/07_plot_methylation_qc.py \
      --input "$input" --output-dir "$results/qc"
    touch "$route/model.COMPLETE"
  done
done
printf 'Recovered all six VMR post-training outputs under %s\n' "$root"
