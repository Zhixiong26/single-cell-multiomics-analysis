#!/usr/bin/env bash

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
output="${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi/vmr_dmr/shared/dmr"
cd "$project"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" \
  Scripts/Methylvi/vmr_dmr/01_prepare_all_unique_pooled_dmrs.py \
  --pairwise-summary "${SCMO_POOLED_DMR_SUMMARY:?set SCMO_POOLED_DMR_SUMMARY}" \
  --blacklist "${SCMO_BLACKLIST:?set SCMO_BLACKLIST}" \
  --chrom-sizes "${SCMO_CHROM_SIZES:?set SCMO_CHROM_SIZES}" \
  --blacklist-fraction "${SCMO_BLACKLIST_FRACTION:?}" --raw-p "${SCMO_HYPO_RAW_P:?}" --min-abs-diff "${SCMO_HYPO_MIN_ABS_DIFF:?}" \
  --sort-threads "${SLURM_CPUS_PER_TASK:-8}" --output-dir "$output"
touch "$output/prepare.COMPLETE"
