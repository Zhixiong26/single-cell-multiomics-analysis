#!/usr/bin/env bash

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
root="${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi/vmr_dmr"
cd "$project"
args=(--root "$root" --expected-cells "${SCMO_EXPECTED_CELLS:?set SCMO_EXPECTED_CELLS}")
for value in ${SCMO_VMR_THRESHOLDS:?set SCMO_VMR_THRESHOLDS}; do args+=(--threshold "$value"); done
for value in ${SCMO_FEATURE_TARGETS:?set SCMO_FEATURE_TARGETS}; do args+=(--feature-target "$value"); done
${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python \
  Scripts/Methylvi/vmr_dmr/03_summarize_vmr_dmr_models.py "${args[@]}"
touch "$root/workflow.COMPLETE"
