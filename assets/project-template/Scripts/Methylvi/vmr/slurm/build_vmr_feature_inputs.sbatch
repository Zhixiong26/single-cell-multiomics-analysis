#!/usr/bin/env bash
# Build one VMR threshold's shared feature input for the largest feature target
# and derive the nested smaller-feature H5MU from it.
# Usage: sbatch [resource options] build_vmr_feature_inputs.sbatch THRESHOLD
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
#SBATCH --job-name=scmo_mvi_vmr_build

set -Eeuo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
if (( $# != 1 )); then echo "Usage: sbatch $0 THRESHOLD" >&2; exit 2; fi
threshold=$1
read -r -a vmr_thresholds <<< "${SCMO_VMR_THRESHOLDS:?set SCMO_VMR_THRESHOLDS}"
[[ " ${vmr_thresholds[*]} " == *" $threshold "* ]] || { echo "Unsupported threshold: $threshold" >&2; exit 3; }
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/vmr
# The route configuration derives its result root from the variance branch, so
# the requested threshold is exported before the configuration is sourced.
export VMR_METHSCAN_VARIANCE="$threshold"
source "$script_dir/00_vmr_methylvi_config.sh"
route_dir="$VMR_RESULTS_ROOT"
export VMR_MVI_INPUT="$route_dir/features_${VMR_TARGET_FEATURES}/input.h5mu"
export VMR_MVI_RESULTS="$route_dir/features_${VMR_TARGET_FEATURES}/results"
export SCMO_MVI_INPUT="$VMR_MVI_INPUT" SCMO_MVI_RESULTS="$VMR_MVI_RESULTS"
export VMR_THREADS="${SLURM_CPUS_PER_TASK:-$VMR_THREADS}"
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" prepare
bash "$script_dir/run.sh" build
python_exe="$SCMO_METHYLVI_ENV/bin/python"
for target in $VMR_FEATURE_TARGETS; do
  [[ "$target" == "$VMR_TARGET_FEATURES" ]] && continue
  "$python_exe" "$script_dir/../shared/06_subset_methylvi_input.py" \
    --input "$VMR_MVI_INPUT" \
    --output "$route_dir/features_${target}/input.h5mu" --features "$target"
done
touch "$route_dir/feature_inputs.COMPLETE"
