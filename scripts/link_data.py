#!/usr/bin/env python3
"""Materialize a project's declared input data under its own Data/ directory.

Every input an analysis reads has to be reachable from the project as `Data/<name>`:
a project whose configuration points at a shared path somewhere else is not
reproducible from its own directory, and it silently breaks the day that path moves.
The `data.sources` block of `config/project.yaml` declares where each input comes
from -- a `path` to link, or a `url` to fetch and checksum.

Planning is the default; it only reads. `--execute` creates the links and performs
the downloads, and it never overwrites an entry that is not the one this source
declares: which copy of the data is real is a human decision, not this tool's.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from _common import (
    WorkflowError, apply_data_sources, load_project, plan_data_sources, stored_data_sources,
)


STATE_NOTES = {
    "ready": "in place and readable from this host",
    "link": "will link to the declared target",
    "dangling": "links to the declared target, which this host cannot see",
    "present": "a file is already here; its sha256 will be verified",
    "download": "will be fetched and checksummed",
    "conflict": "something else is already here; this tool will not overwrite it",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--execute", action="store_true",
                        help="create the links and perform the downloads; without it nothing is written")
    parser.add_argument("--only", action="append", default=[], metavar="NAME",
                        help="restrict the run to these Data/ entries; repeatable")
    args = parser.parse_args()
    root = args.project.resolve()
    entries = stored_data_sources(load_project(root))
    if args.only:
        unknown = sorted(set(args.only) - {entry["name"] for entry in entries})
        if unknown:
            raise WorkflowError("no such data source: %s" % ", ".join(unknown))
        entries = [entry for entry in entries if entry["name"] in set(args.only)]
    if not entries:
        print("No data.sources are declared in config/project.yaml; nothing to link.")
        return 0
    plans = plan_data_sources(root, entries)
    for plan in plans:
        print("%-16s %-9s %s" % (plan["name"], plan["state"], plan["destination"]))
        print("%-16s %-9s %s" % ("", "from", plan["url"] or plan["target"]))
        print("%-16s %-9s %s" % ("", "", STATE_NOTES[plan["state"]]))
    if not args.execute:
        actionable = [plan for plan in plans if plan["state"] not in {"ready"}]
        print("\nPlan only: %d of %d entries need work. Re-run with --execute to apply it."
              % (len(actionable), len(plans)))
        return 0
    for result in apply_data_sources(root, plans, download=True):
        print("%-16s %s" % (result["name"], result["outcome"]))
    print("\nDone. Validate the result with `python tools/validate_project.py --project . --mode full`.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
