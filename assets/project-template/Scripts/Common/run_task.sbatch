#!/usr/bin/env bash
# Resources and dependencies are injected by tools/submit_workflow.py.
set -Eeuo pipefail
if (( $# == 0 )); then
  echo "Usage: $0 COMMAND [ARG ...]" >&2
  exit 2
fi
if [[ -n "${SCMO_RESOURCE_USAGE:-}" ]]; then
  mkdir -p "$(dirname "$SCMO_RESOURCE_USAGE")"
  # GNU time is a separate package and is not installed everywhere. Nothing in
  # the workflow reads the file it writes -- `sacct` supplies the recorded
  # resource usage, and inspect_run.py reads that -- so a task must not fail for
  # want of it. Say so rather than staying silent: a missing
  # resource_usage/<task>.txt is otherwise indistinguishable from a task that
  # never ran.
  if [[ -x /usr/bin/time ]]; then
    exec /usr/bin/time -v -o "$SCMO_RESOURCE_USAGE" "$@"
  fi
  echo "NOTE: /usr/bin/time is absent; running without verbose resource accounting." >&2
  echo "      ${SCMO_RESOURCE_USAGE} will not be written." >&2
fi
exec "$@"
