#!/usr/bin/env python3
"""Execute one planned task command with provenance and completion markers."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from _common import WorkflowError, load_project, signature, write_json  # noqa: E402


def render_command(command, variables):
    tokens = shlex.split(command) if isinstance(command, str) else list(command)
    return [str(token).format(**variables) for token in tokens]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--task", required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    run_dir = root / ".workflow" / "runs" / args.run_id
    plan = json.loads((run_dir / "plan.json").read_text(encoding="utf-8"))
    matches = [item for item in plan["tasks"] if item["id"] == args.task]
    if len(matches) != 1:
        raise WorkflowError("task is not uniquely defined in plan: %s" % args.task)
    item = matches[0]
    cfg = load_project(root)
    commands = cfg["analysis"].get("task_commands") or {}
    command = commands.get(args.task)
    if command is None and args.task == "workflow_summary":
        command = ["/bin/true"]
    if command is None:
        for prefix, value in commands.items():
            if prefix.endswith("*") and args.task.startswith(prefix[:-1]):
                command = value
                break
    if command is None:
        raise WorkflowError(
            "analysis.task_commands does not define %s; generate/adapt the task command before execution" % args.task
        )
    task_dir = run_dir / "tasks" / args.task
    task_dir.mkdir(parents=True, exist_ok=True)
    variables = {
        "project": str(root), "run_id": args.run_id, "run_dir": str(run_dir),
        "task": args.task, "task_dir": str(task_dir), "python": sys.executable,
    }
    variables.update({key: str(value) for key, value in item.get("parameters", {}).items()})
    rendered = render_command(command, variables)
    state = {
        "task": args.task, "status": "running", "command": rendered,
        "started_at": datetime.now(timezone.utc).isoformat(),
        "input_signature": plan["input_signature"],
        "task_signature": signature({"command": rendered, "parameters": item.get("parameters", {})}),
        "allocated_cpus": os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("SCMO_CPUS")),
        "allocated_memory_mb": os.environ.get("SLURM_MEM_PER_NODE", os.environ.get("SCMO_MEMORY_MB")),
    }
    write_json(task_dir / "task_status.json", state)
    code = subprocess.call(rendered, cwd=str(root))
    state["status"] = "complete" if code == 0 else "failed"
    state["return_code"] = code
    state["finished_at"] = datetime.now(timezone.utc).isoformat()
    write_json(task_dir / "task_status.json", state)
    if code == 0:
        (task_dir / "task.COMPLETE").touch()
    return code


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
