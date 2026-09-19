#!/usr/bin/env bash
# Run the complete VMR route for the selected variance branch: input
# preparation, feature build, training, plots, and supervised UMAP.
# Usage: sbatch [resource options] run_methylvi_vmrs.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
# The branch is selected by VMR_METHSCAN_VARIANCE in the route configuration.
#SBATCH --job-name=scmo_mvi_vmr_all

set -euo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/vmr
source "$script_dir/00_vmr_methylvi_config.sh"
export VMR_THREADS="${SLURM_CPUS_PER_TASK:-$VMR_THREADS}"
export VMR_NUMERIC_THREADS=1
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" all
