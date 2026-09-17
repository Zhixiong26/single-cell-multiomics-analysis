#!/usr/bin/env bash
# Usage: sbatch train_vmr_dmr_model.sbatch THRESHOLD 10000|30000

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: $0 0.01|0.02|0.05 10000|30000" >&2; exit 2; fi
case "$1" in 0.01|0.02|0.05) ;; *) exit 3;; esac
case "$2" in 10000|30000) ;; *) exit 4;; esac
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
route="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR/vmr_dmr/var_${1}/features_${2}"
input="$route/input.h5mu"; results="$route/results"
[[ -f "$route/input.COMPLETE" && -s "$input" ]] || exit 5
export SCMO_MVI_INPUT="$input" SCMO_MVI_RESULTS="$results" VMR_MVI_RESULTS="$results"
export SCMO_BATCH_KEY=sample_id SCMO_CELLTYPE_KEY=cell_type SCMO_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export SCMO_SEED=0 SCMO_EPOCHS=500 SCMO_BATCH_SIZE=32
export SCMO_SUPERVISED_TARGET_KEY=cell_type SCMO_SUPERVISED_TARGET_WEIGHTS="0.2 0.5 0.7 0.9"
export SCMO_SUPERVISED_NEIGHBORS=15 SCMO_SUPERVISED_MIN_DIST=0.5 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python_exe=${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python
cd "$project"
"$python_exe" Scripts/Methylvi/shared/04_train_methylvi.py
"$python_exe" Scripts/Methylvi/vmr/03_plot_vmr_methylvi_umap.py \
  --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
"$python_exe" Scripts/Methylvi/shared/05_plot_supervised_umap.py
"$python_exe" Scripts/Methylvi/shared/07_plot_methylation_qc.py --input "$input" --output-dir "$results/qc"
touch "$route/model.COMPLETE"
