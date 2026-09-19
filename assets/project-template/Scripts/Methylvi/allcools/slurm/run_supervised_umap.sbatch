#!/usr/bin/env bash
# Recompute the supervised UMAP for the trained ALLCools model.
# Usage: sbatch [resource options] run_supervised_umap.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
#SBATCH --job-name=scmo_mvi_supumap

set -euo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/allcools
source "$script_dir/00_methylvi_config.sh"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-$SCMO_THREADS}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
bash "$script_dir/run.sh" supervised
