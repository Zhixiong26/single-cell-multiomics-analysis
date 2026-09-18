#!/usr/bin/env bash
# Usage: sbatch combine_vmr_dmr_input.sbatch THRESHOLD 10000|30000

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: $0 0.01|0.02|0.05 10000|30000" >&2; exit 2; fi
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
old="${SCMO_BASE_VMR_INPUT_ROOT:-${SCMO_RESULT_DIR:?}/methylvi/vmr}/var_${1}/features_${2}/input.h5mu"
root="$SCMO_RESULT_DIR/methylvi/vmr_dmr"
dmr="$root/shared/dmr/all_unique_hypo_DMRs.input.h5mu"
out="$root/vmr_dmr/var_${1}/features_${2}/input.h5mu"
[[ -f "$root/shared/dmr/counts.COMPLETE" ]] || exit 5
cd "$project"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" \
  Scripts/Methylvi/vmr_dmr/02_combine_vmr_and_dmr_inputs.py \
  --vmr-input "$old" --dmr-input "$dmr" --output "$out"
touch "$(dirname "$out")/input.COMPLETE"
