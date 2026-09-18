#!/usr/bin/env bash

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
root="${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi/vmr_dmr"
bed="$root/shared/dmr/all_unique_hypo_DMRs.merged.bed"
input="$root/shared/dmr/all_unique_hypo_DMRs.input.h5mu"
work="$root/shared/dmr/count_rows"
allc=${SCMO_ALLC_TABLE:?set SCMO_ALLC_TABLE}
annotation=${SCMO_ANNOTATION:?set SCMO_ANNOTATION}
[[ -f "$root/shared/dmr/prepare.COMPLETE" ]] || { echo "DMR preparation incomplete" >&2; exit 3; }
features=$(wc -l < "$bed")
(( features >= 2 )) || { echo "Too few merged DMRs: $features" >&2; exit 4; }
cd "$project"
VMR_FILTERED_BED="$bed" VMR_ALLC_TABLE="$allc" VMR_ANNOTATION="$annotation" \
VMR_COUNT_ROWS="$work" VMR_MVI_INPUT="$input" VMR_THREADS="${SLURM_CPUS_PER_TASK:-40}" \
VMR_MIN_COVERED_PERCENT=0 VMR_TARGET_FEATURES="$features" \
  bash "$project/Scripts/Methylvi/vmr/run.sh" build
touch "$root/shared/dmr/counts.COMPLETE"
