#!/usr/bin/env bash
# Resources and dependencies are injected by tools/submit_workflow.py.
set -Eeuo pipefail
if (( $# == 0 )); then
  echo "Usage: $0 COMMAND [ARG ...]" >&2
  exit 2
fi
if [[ -n "${SCMO_RESOURCE_USAGE:-}" ]]; then
  mkdir -p "$(dirname "$SCMO_RESOURCE_USAGE")"
  exec /usr/bin/time -v -o "$SCMO_RESOURCE_USAGE" "$@"
fi
exec "$@"
