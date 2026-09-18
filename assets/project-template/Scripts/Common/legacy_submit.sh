#!/usr/bin/env bash
set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
routes=${1:?routes required}
run_id=${2:?usage: legacy wrapper RUN_ID}
python=${SCMO_ORCHESTRATOR_PYTHON:-$(command -v python3)}
plan="$project/.workflow/runs/$run_id/plan.json"
if [[ ! -f "$plan" ]]; then
  "$python" "$project/workflow.py" plan --routes "$routes" --run-id "$run_id"
fi
exec "$python" "$project/workflow.py" submit --run-id "$run_id"
