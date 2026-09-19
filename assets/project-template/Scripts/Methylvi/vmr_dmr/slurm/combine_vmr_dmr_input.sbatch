#!/usr/bin/env bash
# Join one VMR feature input with the pooled-DMR input.
# Usage: sbatch [resource options] combine_vmr_dmr_input.sbatch THRESHOLD FEATURES
# Resources and log paths are injected by the submitter, not declared here.
# The experiment configuration supplies the run root and the model ladder; it
# requires SCMO_MAX_FEATURES.
#SBATCH --job-name=scmo_vmr_dmr_join

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: sbatch $0 THRESHOLD FEATURES" >&2; exit 2; fi
threshold=$1
features=$2
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/vmr_dmr
source "$script_dir/00_vmr_dmr_config.sh"
[[ " ${SCMO_VMR_DMR_THRESHOLDS[*]} " == *" $threshold "* ]] || { echo "Unsupported threshold: $threshold" >&2; exit 3; }
[[ " ${SCMO_VMR_DMR_FEATURE_TARGETS[*]} " == *" $features "* ]] || { echo "Unsupported feature target: $features" >&2; exit 4; }
old=$(vmr_input "$threshold" "$features")
dmr="$SCMO_VMR_DMR_DMR_INPUT"
out="$(leaf_root "$threshold" "$features")/input.h5mu"
[[ -f "$SCMO_VMR_DMR_SHARED/counts.COMPLETE" ]] || { echo "DMR counts are incomplete: $SCMO_VMR_DMR_SHARED" >&2; exit 5; }
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" "$script_dir/02_combine_vmr_and_dmr_inputs.py" \
  --vmr-input "$old" --dmr-input "$dmr" --output "$out"
touch "$(dirname "$out")/input.COMPLETE"
