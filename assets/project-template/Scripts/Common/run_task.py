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
from _common import (WorkflowError, load_project, signature, task_override,
                     validate_recorded_outputs, write_json)  # noqa: E402


def render_command(command, variables):
    tokens = shlex.split(command) if isinstance(command, str) else list(command)
    return [str(token).format(**variables) for token in tokens]


def validate_outputs(item, variables):
    rendered = []
    for value in item.get("outputs", []):
        path = Path(str(value).format(**variables))
        rendered.append(str(path))
    valid, reason = validate_recorded_outputs(rendered)
    if not valid:
        raise WorkflowError(reason or "declared task output evidence is invalid")
    return rendered


def summarize_workflow(plan, run_dir):
    summary_task = next(item for item in plan["tasks"] if item["id"] == "workflow_summary")
    failures = []
    for dependency in summary_task["dependencies"]:
        dependency_plan = next(item for item in plan["tasks"] if item["id"] == dependency)
        task_dir = run_dir / "tasks" / dependency
        marker = task_dir / "task.COMPLETE"
        status_path = task_dir / "task_status.json"
        if not marker.is_file() or not status_path.is_file():
            failures.append("%s lacks completion evidence" % dependency)
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        expected_signature = signature({"command": status.get("command", []),
                                        "parameters": dependency_plan.get("parameters", {})})
        outputs_valid, output_reason = validate_recorded_outputs(status.get("outputs", []))
        if (status.get("status") != "complete" or status.get("input_signature") != plan.get("input_signature")
                or status.get("code_signature") != plan.get("code_signature")
                or status.get("task_signature") != expected_signature
                or not outputs_valid):
            failures.append("%s has invalid status/signature/output evidence: %s" %
                            (dependency, output_reason or "signature mismatch"))
    if failures:
        raise WorkflowError("workflow summary failed: %s" % "; ".join(failures))
    payload = {"status": "complete", "tasks": len(plan["tasks"]), "input_signature": plan["input_signature"]}
    write_json(run_dir / "workflow_summary.json", payload)
    return [str(run_dir / "workflow_summary.json")]


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
    command = task_override(commands, args.task)
    internal_summary = args.task == "workflow_summary" and command is None
    if command is None and not internal_summary:
        command = [sys.executable, str(root / "Scripts" / "Common" / "task_adapter.py"),
                   "--project", str(root), "--run-id", args.run_id, "--task", args.task,
                   "--task-dir", str(run_dir / "tasks" / args.task)]
    task_dir = run_dir / "tasks" / args.task
    task_dir.mkdir(parents=True, exist_ok=True)
    variables = {
        "project": str(root), "run_id": args.run_id, "run_dir": str(run_dir),
        "task": args.task, "task_dir": str(task_dir), "python": sys.executable,
        "result_dir": plan.get("result_dir", str(root / "Results" / "runs" / args.run_id)),
    }
    variables.update({key: str(value) for key, value in item.get("parameters", {}).items()})
    rendered = ["<internal:workflow_summary>"] if internal_summary else render_command(command, variables)
    state = {
        "task": args.task, "status": "running", "command": rendered,
        "started_at": datetime.now(timezone.utc).isoformat(),
        "input_signature": plan["input_signature"], "code_signature": plan.get("code_signature"),
        "task_signature": signature({"command": rendered, "parameters": item.get("parameters", {})}),
        "allocated_cpus": os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("SCMO_CPUS")),
        "allocated_memory_mb": os.environ.get("SLURM_MEM_PER_NODE", os.environ.get("SCMO_MEMORY_MB")),
    }
    write_json(task_dir / "task_status.json", state)
    try:
        if internal_summary:
            outputs = summarize_workflow(plan, run_dir)
            code = 0
        else:
            process_code = subprocess.call(rendered, cwd=str(root))
            state["process_return_code"] = process_code
            code = process_code
            outputs = validate_outputs(item, variables) if process_code == 0 else []
    except (OSError, ValueError, WorkflowError) as exc:
        code, outputs = 1, []
        state["failure_stage"] = "output_validation" if state.get("process_return_code") == 0 else "execution"
        state["error"] = str(exc)
    state["status"] = "complete" if code == 0 else "failed"
    state["return_code"] = code
    state["outputs"] = outputs
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
