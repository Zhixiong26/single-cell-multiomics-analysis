#!/usr/bin/env bash
# Reduced ALLCools-route smoke run: a small cell subset and reduced epochs in
# its own root, so the full run's inputs are never touched.
# Usage: sbatch [resource options] run_methylvi_allcools_smoke.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
# The smoke stage ran verify -> allcools -> build -> train with the overrides
# below; those are the stage names this template's route runner exposes.
#SBATCH --job-name=scmo_mvi_allc_smoke

set -euo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/allcools
source "$script_dir/00_methylvi_config.sh"
smoke_root="${SCMO_SMOKE_ROOT:-${SCMO_MVI_ROOT:?set SCMO_MVI_ROOT}/smoke}"
export SCMO_MVI_ROOT="$smoke_root"
export SCMO_ALLCOOLS_ROOT="$smoke_root/allcools_features"
export SCMO_ALLC_DIR="$SCMO_ALLCOOLS_ROOT/input_allc"
export SCMO_ALLC_TABLE="$SCMO_ALLCOOLS_ROOT/selected_cells.allc.tsv"
export SCMO_MCDS="$SCMO_ALLCOOLS_ROOT/mcg.mcds"
export SCMO_ALLCOOLS_H5AD="$SCMO_ALLCOOLS_ROOT/clustered.h5ad"
export SCMO_MVI_INPUT="$smoke_root/methylvi_input.h5mu"
export SCMO_MVI_RESULTS="$smoke_root/results"
export SCMO_MAX_CELLS="${SCMO_SMOKE_CELLS:-100}"
export SCMO_BALANCED_COHORTS=1
export SCMO_EPOCHS="${SCMO_SMOKE_EPOCHS:-2}"
export SCMO_NUMERIC_THREADS=1
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-$SCMO_THREADS}"
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" verify
bash "$script_dir/run.sh" allcools
bash "$script_dir/run.sh" build
bash "$script_dir/run.sh" train
