#!/usr/bin/env python3
"""Inspect planned/submitted task state and emit a machine-readable run summary."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from _common import (
    TERMINAL_BAD_STATES, WorkflowError, refresh_run_log, signature,
    validate_recorded_outputs, write_json,
)


TERMINAL_OK = {"COMPLETED", "complete"}
# Aliased from the shared set: the run's status verdict and the "failed" count in
# the Report record are computed from the same states, so they cannot disagree.
TERMINAL_BAD = TERMINAL_BAD_STATES


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
    submission_by_task = {item["task"]: item for item in submissions}
    task_rows = []
    for planned in plan.get("tasks", []):
        item = submission_by_task.get(planned["id"])
        job_id = str(item.get("job_id", "")) if item else ""
        scheduler_state = usage.get(job_id, {}).get(
            "state", item.get("status", "unknown") if item else "not_submitted"
        )
        task_dir = run_dir / "tasks" / planned["id"]
        marker = task_dir / "task.COMPLETE"
        status_path = task_dir / "task_status.json"
        evidence_valid, reason = False, None
        if marker.is_file() and status_path.is_file():
            try:
                status_record = json.loads(status_path.read_text(encoding="utf-8"))
                expected_signature = signature({
                    "command": status_record.get("command", []),
                    "parameters": planned.get("parameters", {}),
                })
                outputs_valid, output_reason = validate_recorded_outputs(status_record.get("outputs", []))
                evidence_valid = (
                    status_record.get("status") == "complete"
                    and status_record.get("return_code") == 0
                    and status_record.get("input_signature") == plan.get("input_signature")
                    and status_record.get("code_signature") == plan.get("code_signature")
                    and status_record.get("task_signature") == expected_signature
                    and outputs_valid
                )
                if not evidence_valid:
                    reason = output_reason or "completion evidence or signature mismatch"
            except (OSError, ValueError, TypeError):
                reason = "invalid task_status.json"
        elif marker.is_file() or status_path.is_file():
            reason = "incomplete marker/status evidence"
        state = "complete" if evidence_valid else scheduler_state
        task_rows.append({
            "task": planned["id"], "job_id": job_id, "state": state,
            "complete_marker": marker.is_file(), "evidence_valid": evidence_valid,
            "evidence_error": reason, "resource_usage": usage.get(job_id),
        })
    if any(row["state"] in TERMINAL_BAD for row in task_rows):
        status = "failed"
    elif task_rows and all(row["evidence_valid"] for row in task_rows):
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
    else:
        marker = run_dir / "workflow.COMPLETE"
        if marker.exists():
            marker.unlink()
    # Deliberately after both of the above: the record points at run_summary.json
    # and the refresh rescans every run, so this run's summary has to be on disk
    # first, and the document must never claim more than the machine state does.
    _refresh_docs(Path(project).resolve())
    return result


def _refresh_docs(project: Path) -> None:
    """Record this run in the root Report.md, best effort.

    The exit code encodes the state of the run, not of a derived document: a
    read-only checkout or a full disk must not turn a validated run into a
    reported failure. update_report.py is the loud path, and it does exit
    non-zero when the document cannot be written, so a failure is never silent
    to someone who wants it to be fatal.
    """
    try:
        outcome = refresh_run_log(project)
    except (WorkflowError, OSError, ValueError) as exc:
        print(
            "WARNING: could not update %s: %s\n"
            "         recover with: python tools/update_report.py --project %s"
            % (project / "Report.md", exc, project),
            file=sys.stderr,
        )
        return
    for warning in outcome.get("warnings") or []:
        print("WARNING: run summary not readable: %s" % warning, file=sys.stderr)


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
