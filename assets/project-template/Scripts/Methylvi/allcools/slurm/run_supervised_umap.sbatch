#!/usr/bin/env bash
#SBATCH --ntasks=1

set -euo pipefail
export SCMO_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd ${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
bash Scripts/Methylvi/allcools/run.sh supervised
