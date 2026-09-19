#!/usr/bin/env bash
# Emit sbatch resource flags for one scheduler profile.
#
# Stage scripts carry no site-specific #SBATCH directives. Resources are read
# from config/scheduler.yaml and passed on the sbatch command line, which is the
# same contract tools/submit_workflow.py uses for the managed run DAG.
#
# Usage: scheduler_flags.sh PROFILE
# Prints: --partition P --cpus-per-task N --mem M --time T [--account A]
set -Eeuo pipefail
# SCMO_PROJECT_ROOT wins when the caller knows it; otherwise locate the project
# from this script's own path (Scripts/Common/ -> project root). That keeps the
# by-hand call used for exploring a stage working without any setup.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project=${SCMO_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
profile=${1:?usage: scheduler_flags.sh PROFILE}
python=${SCMO_ORCHESTRATOR_PYTHON:-$(command -v python3)}

"$python" - "$project" "$profile" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
profile_name = sys.argv[2]
sys.path.insert(0, str(root / "tools"))
from _common import load_project  # noqa: E402

scheduler = load_project(root).get("scheduler") or {}
profiles = scheduler.get("profiles") or {}
profile = profiles.get(profile_name) or profiles.get("default")
if not profile:
    raise SystemExit("ERROR: no scheduler profile named %r and no default profile" % profile_name)
target = profile.get("target") or {}
missing = [key for key in ("cpus", "memory") if not target.get(key)]
if missing:
    raise SystemExit("ERROR: profile %r is missing target.%s" % (profile_name, ", target.".join(missing)))

flags = ["--cpus-per-task", str(target["cpus"]), "--mem", str(target["memory"])]
if target.get("time"):
    flags += ["--time", str(target["time"])]
# A profile may name a partition; otherwise let the cluster default apply.
if profile.get("partition"):
    flags += ["--partition", str(profile["partition"])]
if scheduler.get("account"):
    flags += ["--account", str(scheduler["account"])]
print(" ".join(flags))
PY
