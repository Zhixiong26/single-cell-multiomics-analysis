#!/usr/bin/env python3
"""Inspect planned/submitted task state and emit a machine-readable run summary."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from _common import WorkflowError, write_json


TERMINAL_OK = {"COMPLETED", "complete"}
TERMINAL_BAD = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "failed"}


def sacct_usage(job_ids: list[str]) -> dict[str, dict]:
    numeric = [value for value in job_ids if value.isdigit()]
    if not numeric:
        return {}
    try:
        output = subprocess.check_output([
            "sacct", "-n", "-P", "-j", ",".join(numeric),
            "-o", "JobIDRaw,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,TotalCPU,ExitCode"
        ], stderr=subprocess.STDOUT).decode("utf-8", "replace")
    except (OSError, subprocess.SubprocessError):
        return {}
    result = {}
    for line in output.splitlines():
        fields = line.split("|")
        if len(fields) >= 2 and "." not in fields[0]:
            result[fields[0]] = {
                "state": fields[1].split()[0].split("+")[0],
                "elapsed": fields[2] if len(fields) > 2 else "",
                "allocated_cpus": fields[3] if len(fields) > 3 else "",
                "requested_memory": fields[4] if len(fields) > 4 else "",
                "max_rss": fields[5] if len(fields) > 5 else "",
                "total_cpu": fields[6] if len(fields) > 6 else "",
                "exit_code": fields[7] if len(fields) > 7 else "",
            }
    return result


def inspect(project: Path, run_id: str) -> dict:
    run_dir = Path(project).resolve() / ".workflow" / "runs" / run_id
    plan_path = run_dir / "plan.json"
    if not plan_path.is_file():
        raise WorkflowError("run plan is absent: %s" % plan_path)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    submissions_path = run_dir / "submissions.json"
    submissions = json.loads(submissions_path.read_text(encoding="utf-8")) if submissions_path.is_file() else []
    usage = sacct_usage([str(item.get("job_id", "")) for item in submissions])
    write_json(run_dir / "resource_usage" / "sacct.json", usage)
    task_rows = []
    for item in submissions:
        job_id = str(item.get("job_id", ""))
        state = usage.get(job_id, {}).get("state", item.get("status", "unknown"))
        marker = run_dir / "tasks" / item["task"] / "task.COMPLETE"
        if marker.is_file():
            state = "complete"
        task_rows.append({
            "task": item["task"], "job_id": job_id, "state": state,
            "complete_marker": marker.is_file(), "resource_usage": usage.get(job_id),
        })
    if any(row["state"] in TERMINAL_BAD for row in task_rows):
        status = "failed"
    elif task_rows and all(row["state"] in TERMINAL_OK for row in task_rows):
        status = "complete"
    elif submissions:
        status = "in_progress"
    else:
        status = "planned"
    result = {
        "schema_version": 1, "run_id": run_id, "status": status,
        "checked_at": datetime.now(timezone.utc).isoformat(), "tasks": task_rows,
        "input_signature": plan.get("input_signature"),
    }
    write_json(run_dir / "run_summary.json", result)
    if status == "complete":
        (run_dir / "workflow.COMPLETE").touch()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    result = inspect(args.project, args.run_id)
    if args.json_out:
        write_json(args.json_out, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
