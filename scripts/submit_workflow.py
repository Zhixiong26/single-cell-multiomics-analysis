#!/usr/bin/env python3
"""Submit or locally execute a planned DAG with a fresh resource check per task."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from _common import (HOST_ROLE_SUBMIT, LOGIN_EXECUTION_VARIABLE, WorkflowError, available_partitions,
                     effective_backend, host_role, host_role_source, load_project, login_execution_ack,
                     task_is_implemented, validate_recorded_outputs, write_json)
from inspect_resources import inspect
from validate_project import validate


def load_plan(project: Path, run_id: str) -> tuple[Path, dict]:
    run_dir = Path(project).resolve() / ".workflow" / "runs" / run_id
    path = run_dir / "plan.json"
    if not path.is_file():
        raise WorkflowError("plan does not exist: %s" % path)
    return run_dir, json.loads(path.read_text(encoding="utf-8"))


def environment_python(project: Path, stage: str = "orchestrator") -> str:
    path = Path(project).resolve() / "config" / "environments.tsv"
    import csv
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    preferred = [stage]
    if stage == "orchestrator":
        preferred.append("common")
    for wanted in preferred:
        for row in rows:
            if row["stage"] == wanted and row["python"].strip():
                return row["python"].strip()
    for row in rows:
        if row["stage"] in {"orchestrator", "common"} and row["python"].strip():
            return row["python"].strip()
    raise WorkflowError("no Python interpreter is configured for stage %s" % stage)


def task_command(project: Path, run_id: str, item: dict) -> list[str]:
    root = Path(project).resolve()
    return [
        environment_python(root, item.get("environment", "orchestrator")),
        str(root / "Scripts" / "Common" / "run_task.py"),
        "--project", str(root), "--run-id", run_id, "--task", item["id"],
    ]


def completed_evidence(run_dir: Path, task_id: str, input_signature: str) -> bool:
    task_dir = run_dir / "tasks" / task_id
    marker, status_path = task_dir / "task.COMPLETE", task_dir / "task_status.json"
    if not marker.is_file() or not status_path.is_file():
        return False
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    outputs_valid, _ = validate_recorded_outputs(status.get("outputs", []))
    return (status.get("status") == "complete" and status.get("return_code") == 0
            and status.get("input_signature") == input_signature and outputs_valid)


def slurm_job_state(job_id: str) -> str:
    """Return the current Slurm state or fail closed when it cannot be resolved."""
    try:
        queued = subprocess.check_output(
            ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
            stderr=subprocess.STDOUT,
        ).decode("utf-8").strip().splitlines()
    except (OSError, subprocess.SubprocessError) as exc:
        raise WorkflowError("cannot query squeue for job %s: %s" % (job_id, exc))
    if queued:
        return queued[0].strip().upper().split("+", 1)[0]
    try:
        history = subprocess.check_output(
            ["sacct", "-n", "-X", "-j", str(job_id), "--format=State", "--parsable2"],
            stderr=subprocess.STDOUT,
        ).decode("utf-8").strip().splitlines()
    except (OSError, subprocess.SubprocessError) as exc:
        raise WorkflowError("cannot query sacct for job %s: %s" % (job_id, exc))
    states = [line.split("|", 1)[0].strip().upper().split("+", 1)[0]
              for line in history if line.split("|", 1)[0].strip()]
    if not states:
        raise WorkflowError("Slurm job %s is absent from both squeue and sacct; refusing duplicate submission" % job_id)
    return states[0]


def submit(project: Path, run_id: str, dry_run: bool = False,
           allow_login_execution: bool = False) -> dict:
    root = Path(project).resolve()
    run_dir, plan = load_plan(root, run_id)
    required_stages = {stage for item in plan.get("tasks", [])
                       for stage in item.get("required_environments", [item.get("environment", "orchestrator")])}
    selected_routes = {name for name, enabled in plan.get("routes", {}).items() if enabled}
    current = validate(root, required_stages=required_stages, selected_routes=selected_routes)
    if current["status"] != "valid":
        raise WorkflowError("pre-submit preflight failed: %s" % "; ".join(current["errors"]))
    if current["input_signature"] != plan["input_signature"]:
        raise WorkflowError("input signature changed after planning; create a new run-id")
    cfg = load_project(root)
    scheduler = cfg["scheduler"]
    role = host_role()
    # The acknowledgement is passed as a value into the environment the task will
    # see, never written into os.environ: an inherited variable would ride along
    # into an sbatch job (--export ALL) and into later runs, and a bypass that
    # propagates itself is not bound to the run it was granted for.
    child_env = os.environ.copy()
    if allow_login_execution:
        child_env[LOGIN_EXECUTION_VARIABLE] = run_id
    backend, backend_source = effective_backend(scheduler, run_id, role=role, environ=child_env)
    if backend == "slurm" and not scheduler.get("partitions"):
        if backend_source == "promoted-from-local":
            observed = available_partitions()
            raise WorkflowError(
                "scheduler.backend is local but this host is a Slurm submit host, so run %s was "
                "promoted to Slurm rather than executed on the login node; that needs a partition "
                "allow-list and scheduler.partitions is empty. List the partitions this project may "
                "use%s." % (run_id, " -- this cluster offers %s" % ", ".join(observed) if observed else ""))
        raise WorkflowError("Slurm backend requires a non-empty scheduler.partitions allow-list")
    if not dry_run:
        validation_path = root / ".workflow" / "validations" / (current["input_signature"] + ".full.json")
        full = None
        if validation_path.is_file():
            try:
                candidate = json.loads(validation_path.read_text(encoding="utf-8"))
                if (candidate.get("status") == "valid" and candidate.get("validation_mode") == "full"
                        and candidate.get("input_signature") == current["input_signature"]):
                    full = candidate
            except (OSError, ValueError):
                pass
        if full is None:
            full = validate(root, mode="full", workers=max(1, int(scheduler.get("validation_workers", 4))),
                            required_stages=required_stages, selected_routes=selected_routes)
            write_json(validation_path, full)
        if full["status"] != "valid" or full["input_signature"] != current["input_signature"]:
            raise WorkflowError("matching full input validation is required before production submission")
    if not dry_run:
        commands = cfg["analysis"].get("task_commands") or {}
        missing = [item["id"] for item in plan["tasks"] if not task_is_implemented(item["id"], commands)]
        if missing:
            raise WorkflowError(
                "production submission is blocked because task commands are not adapted: %s"
                % ", ".join(missing)
            )
    wrapper = root / "Scripts" / "Common" / "run_task.sbatch"
    if backend == "slurm" and not wrapper.is_file():
        raise WorkflowError("Slurm wrapper is absent: %s" % wrapper)

    submissions_path = run_dir / "submissions.json"
    previous = json.loads(submissions_path.read_text(encoding="utf-8")) if submissions_path.is_file() else []
    records = {item["task"]: item for item in previous}
    job_ids = {}
    max_parallel = max(1, int(scheduler.get("max_parallel", 2)))
    limited_profiles = set(scheduler.get("limited_profiles") or ["methscan_branch", "dmr", "feature_builder", "trainer"])
    throttle_slots = [None] * max_parallel
    throttle_index = 0
    for item in plan["tasks"]:
        task_id = item["id"]
        if completed_evidence(run_dir, task_id, plan["input_signature"]):
            records[task_id] = dict(records.get(task_id, {}), task=task_id, status="complete", resumed=True)
            write_json(submissions_path, [records[x["id"]] for x in plan["tasks"] if x["id"] in records])
            continue
        existing = records.get(task_id)
        if existing and existing.get("status") in {"submitted", "running"} and existing.get("job_id"):
            if backend != "slurm":
                raise WorkflowError("non-Slurm task %s has a stale in-progress submission record" % task_id)
            state = slurm_job_state(str(existing["job_id"]))
            existing["scheduler_state"] = state
            if state in {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED", "RESIZING"}:
                job_ids[task_id] = str(existing["job_id"])
                write_json(submissions_path, [records[x["id"]] for x in plan["tasks"] if x["id"] in records])
                continue
            existing["status"] = "retryable"
            existing["previous_job_id"] = str(existing["job_id"])
        snapshot_path = run_dir / "resource_snapshots" / (task_id + ".json")
        snapshot = inspect(root, item["profile"], run_id=run_id, resolved_backend=backend)
        write_json(snapshot_path, snapshot)
        recommendation = snapshot["recommendation"]
        command = task_command(root, run_id, item)
        unresolved = [name for name in item["dependencies"]
                      if name not in job_ids and not completed_evidence(run_dir, name, plan["input_signature"])]
        if unresolved:
            raise WorkflowError("task %s has unresolved dependencies: %s" % (task_id, ", ".join(unresolved)))
        dependencies = [job_ids[name] for name in item["dependencies"] if name in job_ids]
        throttle_dependencies = []
        slot = None
        if backend == "slurm" and item["profile"] in limited_profiles:
            slot = throttle_index % max_parallel
            prior = throttle_slots[slot]
            if prior and prior not in dependencies:
                dependencies.append(prior)
                throttle_dependencies.append(prior)
            throttle_index += 1
        record = {
            "task": task_id, "profile": item["profile"], "dependencies": item["dependencies"],
            "throttle_dependencies": throttle_dependencies,
            "resource_snapshot": str(snapshot_path), "recommendation": recommendation,
            "command": command, "submitted_at": datetime.now(timezone.utc).isoformat(),
            # Which host decided, and which backend really ran. A record that does
            # not say both cannot show that work stayed off the login node, and
            # `backend_declared` differing from `backend_effective` is what a
            # promotion looks like after the fact.
            "host_role": role, "host_role_source": host_role_source(),
            "backend_declared": str(scheduler.get("backend") or "local").strip().lower(),
            "backend_effective": backend, "backend_source": backend_source,
            "login_execution": backend == "local" and role == HOST_ROLE_SUBMIT,
            "login_execution_ack": login_execution_ack(child_env) if backend == "local" else "",
        }
        if dry_run:
            record["status"] = "dry_run"
            record["job_id"] = "dry_%s" % task_id
            job_ids[task_id] = record["job_id"]
        elif backend == "local":
            env = dict(child_env)
            env["SCMO_CPUS"] = str(recommendation["cpus"])
            env["SCMO_MEMORY_MB"] = str(recommendation["memory_mb"])
            code = subprocess.call(command, cwd=str(root), env=env)
            if code:
                record["status"] = "failed"
                record["return_code"] = code
                records[task_id] = record
                write_json(submissions_path, [records[x["id"]] for x in plan["tasks"] if x["id"] in records])
                raise WorkflowError("local task failed: %s" % task_id)
            record["status"] = "complete"
            record["job_id"] = "local_%s" % task_id
            job_ids[task_id] = record["job_id"]
        elif backend == "slurm":
            output = run_dir / "logs" / (task_id + "_%j.out")
            error = run_dir / "logs" / (task_id + "_%j.err")
            output.parent.mkdir(parents=True, exist_ok=True)
            sbatch = [
                "sbatch", "--parsable", "--job-name", "scmo_%s" % task_id[:60],
                "--partition", recommendation["partition"],
                "--cpus-per-task", str(recommendation["cpus"]),
                "--mem", recommendation["memory"], "--time", recommendation["time"],
                "--output", str(output), "--error", str(error),
                "--export", "ALL,SCMO_RESOURCE_USAGE=%s" % (run_dir / "resource_usage" / (task_id + ".txt")),
            ]
            account = scheduler.get("account")
            if account:
                sbatch.extend(["--account", str(account)])
            if dependencies:
                sbatch.extend(["--dependency", "afterok:%s" % ":".join(dependencies)])
            if recommendation.get("pin_node"):
                sbatch.extend(["--nodelist", recommendation["reference_node"]])
            sbatch.extend([str(wrapper)] + command)
            output_text = subprocess.check_output(sbatch, cwd=str(root), stderr=subprocess.STDOUT).decode("utf-8").strip()
            job_id = output_text.split(";", 1)[0]
            record.update({"status": "submitted", "job_id": job_id, "sbatch": sbatch})
            job_ids[task_id] = job_id
            if slot is not None:
                throttle_slots[slot] = job_id
        else:
            raise WorkflowError("unsupported backend: %s" % backend)
        records[task_id] = record
        write_json(submissions_path, [records[x["id"]] for x in plan["tasks"] if x["id"] in records])
    plan["status"] = "dry_run" if dry_run else ("submitted" if backend == "slurm" else "complete")
    plan["job_ids"] = job_ids
    plan["host_role"] = role
    plan["backend_declared"] = str(scheduler.get("backend") or "local").strip().lower()
    plan["backend_effective"] = backend
    plan["backend_source"] = backend_source
    if backend == "local" and role == HOST_ROLE_SUBMIT:
        plan["login_execution_ack"] = login_execution_ack(child_env)
    write_json(run_dir / "plan.json", plan)
    return {"status": plan["status"], "run_dir": str(run_dir), "jobs": job_ids,
            "host_role": role, "backend_effective": backend, "backend_source": backend_source}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-login-execution", action="store_true",
        help="Execute this run's tasks on this host even when it is a Slurm submit host. The "
             "acknowledgement names the run, is recorded in submissions.json and in each task's "
             "status file, and is never inherited by a submitted job or a later run. Use it for "
             "small debugging runs only; prefer running inside an allocation: salloc, or "
             "srun --pty bash.")
    args = parser.parse_args()
    result = submit(args.project, args.run_id, args.dry_run, args.allow_login_execution)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
