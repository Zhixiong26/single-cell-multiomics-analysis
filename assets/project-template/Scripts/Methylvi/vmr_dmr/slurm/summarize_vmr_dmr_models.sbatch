#!/usr/bin/env bash

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
root="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR"
cd "$project"
${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python \
  Scripts/Methylvi/vmr_dmr/03_summarize_vmr_dmr_models.py --root "$root"
touch "$root/workflow.COMPLETE"
