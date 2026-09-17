#!/usr/bin/env python3
"""Convenience dispatcher for generated project workflow tools."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


TOOLS = {
    "validate": "validate_project.py", "plan": "plan_workflow.py",
    "resources": "inspect_resources.py", "submit": "submit_workflow.py",
    "inspect": "inspect_run.py", "report": "update_report.py",
}


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] not in TOOLS:
        print("Usage: workflow.py {%s} [arguments]" % "|".join(sorted(TOOLS)), file=sys.stderr)
        return 2
    root = Path(__file__).resolve().parent
    command = [sys.executable, str(root / "tools" / TOOLS[sys.argv[1]])] + sys.argv[2:]
    if "--project" not in command:
        command.extend(["--project", str(root)])
    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main())
