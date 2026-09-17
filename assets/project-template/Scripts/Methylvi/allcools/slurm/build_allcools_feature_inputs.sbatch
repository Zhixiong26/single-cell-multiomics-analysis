#!/usr/bin/env bash
# Build the shared 30k ALLCools input and derive its nested 10k H5MU.
#SBATCH --ntasks=1

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
source "$project/Scripts/Methylvi/00_experiment_config.sh"
route="$MVI_EXPERIMENT_ROOT/allcools"
shared="$route/shared"
export SCMO_MVI_ROOT="$shared"
export SCMO_ALLCOOLS_ROOT="$shared/allcools_5kb"
export SCMO_ALLC_DIR="$shared/input_allc"
export SCMO_ALLC_TABLE="$shared/selected_cells.allc.tsv"
export SCMO_MCDS="$shared/mcg_5kb.mcds"
export SCMO_ALLCOOLS_H5AD="$route/features_30000/allcools_30k.h5ad"
export SCMO_MVI_INPUT="$route/features_30000/input.h5mu"
export SCMO_MVI_RESULTS="$route/features_30000/results"
export SCMO_TARGET_FEATURES=30000
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-45}"
cd "$project"
bash Scripts/Methylvi/allcools/run.sh prepare
bash Scripts/Methylvi/allcools/run.sh cluster
bash Scripts/Methylvi/allcools/run.sh build
"$MVI_METHYLVI_ENV/bin/python" Scripts/Methylvi/shared/06_subset_methylvi_input.py \
  --input "$route/features_30000/input.h5mu" \
  --output "$route/features_10000/input.h5mu" --features 10000
touch "$route/feature_inputs.COMPLETE"
