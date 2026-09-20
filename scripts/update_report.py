#!/usr/bin/env python3
"""Rebuild the run log in the generated root Report.md from machine-readable run evidence.

The log is a canonical re-render, so this is the same code path the automatic
refresh inside inspect_run.py uses. Reach for it after pruning or archiving a run
directory, after migrating a project written by an older release, or when the
region has been edited by hand and no longer matches the evidence on disk.

With no --run-id the entire region is rebuilt from every
`.workflow/runs/*/run_summary.json`. With --run-id the given run must exist, and
its record is refreshed from its summary; other runs keep their records.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from _common import WorkflowError, refresh_run_log, refresh_stage_run_logs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", help="require this run to exist before rebuilding")
    args = parser.parse_args()
    root = args.project.resolve()
    if args.run_id:
        # A missing summary is an error here rather than a warning: the caller
        # named this run, so silently rebuilding without it would hide the
        # mistake. It is still the same renderer below.
        summary_path = root / ".workflow" / "runs" / args.run_id / "run_summary.json"
        if not summary_path.is_file():
            raise WorkflowError(
                "no run summary at %s; run inspect_run.py --run-id %s first"
                % (summary_path, args.run_id)
            )
    outcome = refresh_run_log(root)
    for warning in outcome.get("warnings") or []:
        print("WARNING: run summary not readable: %s" % warning, file=sys.stderr)
    print("%s (%d record(s), %s)"
          % (outcome["path"], outcome["records"], "rewritten" if outcome["written"] else "unchanged"))
    # The stage logs are the same summaries filtered per stage; a project that
    # has no stage documents yet simply reports none and this prints nothing.
    stage_outcome = refresh_stage_run_logs(root)
    for warning in stage_outcome.get("warnings") or []:
        print("WARNING: run summary not readable: %s" % warning, file=sys.stderr)
    for stage in stage_outcome.get("stages") or []:
        print("%s (%d record(s), %s)"
              % (stage["path"], stage["records"], "rewritten" if stage["written"] else "unchanged"))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
