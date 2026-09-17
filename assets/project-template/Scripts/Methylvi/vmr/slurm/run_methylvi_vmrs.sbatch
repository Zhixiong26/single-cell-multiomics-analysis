#!/usr/bin/env bash
#SBATCH --ntasks=1

set -euo pipefail
export VMR_THREADS="$SLURM_CPUS_PER_TASK"
export VMR_NUMERIC_THREADS=1
cd ${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
bash Scripts/Methylvi/vmr/run.sh all
