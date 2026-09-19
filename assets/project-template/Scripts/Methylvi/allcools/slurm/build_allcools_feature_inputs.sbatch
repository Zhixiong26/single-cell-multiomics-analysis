#!/usr/bin/env bash
# Build the shared ALLCools feature input for the largest feature target and
# derive the nested smaller-feature H5MU from it.
# Usage: sbatch [resource options] build_allcools_feature_inputs.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
#SBATCH --job-name=scmo_mvi_allc_build

set -Eeuo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/allcools
source "$script_dir/00_methylvi_config.sh"
# The experiment configuration requires an explicit maximum feature target; for
# the ALLCools route the largest configured target is that maximum.
export SCMO_MAX_FEATURES="${SCMO_MAX_FEATURES:-$SCMO_TARGET_FEATURES}"
source "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/00_experiment_config.sh"
route_dir="$SCMO_MVI_EXPERIMENT_ROOT/allcools"
shared="$route_dir/shared"
export SCMO_MVI_ROOT="$shared"
export SCMO_ALLCOOLS_ROOT="$shared/allcools_features"
export SCMO_ALLC_DIR="$shared/input_allc"
export SCMO_ALLC_TABLE="$shared/selected_cells.allc.tsv"
export SCMO_MCDS="$shared/mcg.mcds"
export SCMO_ALLCOOLS_H5AD="$route_dir/features_${SCMO_TARGET_FEATURES}/allcools.h5ad"
export SCMO_MVI_INPUT="$route_dir/features_${SCMO_TARGET_FEATURES}/input.h5mu"
export SCMO_MVI_RESULTS="$route_dir/features_${SCMO_TARGET_FEATURES}/results"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-$SCMO_THREADS}"
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" prepare
bash "$script_dir/run.sh" cluster
bash "$script_dir/run.sh" build
python_exe="$SCMO_METHYLVI_ENV/bin/python"
for target in $SCMO_FEATURE_TARGETS; do
  [[ "$target" == "$SCMO_TARGET_FEATURES" ]] && continue
  "$python_exe" "$script_dir/../shared/06_subset_methylvi_input.py" \
    --input "$SCMO_MVI_INPUT" \
    --output "$route_dir/features_${target}/input.h5mu" --features "$target"
done
touch "$route_dir/feature_inputs.COMPLETE"
