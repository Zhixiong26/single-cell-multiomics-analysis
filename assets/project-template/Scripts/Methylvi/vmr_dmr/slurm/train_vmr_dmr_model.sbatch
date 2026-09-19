#!/usr/bin/env bash
# Train and plot one VMR+pooled-DMR model.
# Usage: sbatch [resource options] train_vmr_dmr_model.sbatch THRESHOLD FEATURES
# Resources and log paths are injected by the submitter, not declared here.
# The experiment configuration supplies the run root and the model ladder; it
# requires SCMO_MAX_FEATURES.
#SBATCH --job-name=scmo_vmr_dmr_train

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: sbatch $0 THRESHOLD FEATURES" >&2; exit 2; fi
threshold=$1
features=$2
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi
source "$script_dir/vmr_dmr/00_vmr_dmr_config.sh"
[[ " ${SCMO_VMR_DMR_THRESHOLDS[*]} " == *" $threshold "* ]] || { echo "Unsupported threshold: $threshold" >&2; exit 3; }
[[ " ${SCMO_VMR_DMR_FEATURE_TARGETS[*]} " == *" $features "* ]] || { echo "Unsupported feature target: $features" >&2; exit 4; }
route=$(leaf_root "$threshold" "$features")
input="$route/input.h5mu"; results="$route/results"
[[ -f "$route/input.COMPLETE" && -s "$input" ]] || { echo "Input is not complete: $input" >&2; exit 5; }
export SCMO_MVI_INPUT="$input" SCMO_MVI_RESULTS="$results"
# The VMR plotting entry point reads its default output location from this
# route-specific variable; keep it aligned with the shared MethylVI results.
export VMR_MVI_RESULTS="$results"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-${SCMO_THREADS:-16}}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python_exe="${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python"
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
"$python_exe" "$script_dir/shared/04_train_methylvi.py"
"$python_exe" "$script_dir/vmr/03_plot_vmr_methylvi_umap.py" \
  --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
"$python_exe" "$script_dir/shared/05_plot_supervised_umap.py"
"$python_exe" "$script_dir/shared/07_plot_methylation_qc.py" --input "$input" --output-dir "$results/qc"
touch "$route/model.COMPLETE"
