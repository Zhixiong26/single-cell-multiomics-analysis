#!/usr/bin/env bash

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
output="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR/shared/dmr"
cd "$project"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" \
  Scripts/Methylvi/vmr_dmr/01_prepare_all_unique_pooled_dmrs.py \
  --pairwise-summary "${SCMO_POOLED_DMR_SUMMARY:?set SCMO_POOLED_DMR_SUMMARY}" \
  --blacklist "${SCMO_BLACKLIST:?set SCMO_BLACKLIST}" \
  --blacklist-fraction 0.2 --raw-p 0.01 --min-abs-diff 0.25 \
  --sort-threads "${SLURM_CPUS_PER_TASK:-8}" --output-dir "$output"
touch "$output/prepare.COMPLETE"
