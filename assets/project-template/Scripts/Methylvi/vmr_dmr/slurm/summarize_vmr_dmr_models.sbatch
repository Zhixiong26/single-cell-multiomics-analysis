#!/usr/bin/env bash
# Validate and summarize every VMR+pooled-DMR model of this run.
# Usage: sbatch [resource options] summarize_vmr_dmr_models.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The experiment configuration supplies the run root and the model ladder; it
# requires SCMO_MAX_FEATURES, and the expected cell count comes from the caller.
#SBATCH --job-name=scmo_vmr_dmr_summary

set -Eeuo pipefail
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/vmr_dmr
source "$script_dir/00_vmr_dmr_config.sh"
summarize_args=(--root "$SCMO_VMR_DMR_ROOT" --expected-cells "${SCMO_EXPECTED_CELLS:?set SCMO_EXPECTED_CELLS}")
for threshold in "${SCMO_VMR_DMR_THRESHOLDS[@]}"; do
  summarize_args+=(--threshold "$threshold")
done
for target in "${SCMO_VMR_DMR_FEATURE_TARGETS[@]}"; do
  summarize_args+=(--feature-target "$target")
done
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" "$script_dir/03_summarize_vmr_dmr_models.py" "${summarize_args[@]}"
touch "$SCMO_VMR_DMR_ROOT/workflow.COMPLETE"
