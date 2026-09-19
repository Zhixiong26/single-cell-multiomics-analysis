#!/usr/bin/env bash
# Plot the ALLCools clustering UMAP before MethylVI training.
# Usage: sbatch [resource options] run_allcools_umap_plots.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
#SBATCH --job-name=scmo_mvi_allc_umap

set -euo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/allcools
source "$script_dir/00_methylvi_config.sh"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-$SCMO_THREADS}"
export MPLBACKEND=Agg
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" plots-before
