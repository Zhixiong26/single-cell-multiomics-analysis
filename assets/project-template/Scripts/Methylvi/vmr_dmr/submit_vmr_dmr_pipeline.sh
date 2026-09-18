#!/usr/bin/env bash
set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
exec bash "$project/Scripts/Common/legacy_submit.sh" methylvi_vmr_dmr "${1:?usage: $0 RUN_ID}"
