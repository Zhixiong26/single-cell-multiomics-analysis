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

from _common import WorkflowError, load_project, write_json
from inspect_resources import inspect
from validate_project import validate


def load_plan(project: Path, run_id: str) -> tuple[Path, dict]:
    run_dir = Path(project).resolve() / ".workflow" / "runs" / run_id
    path = run_dir / "plan.json"
    if not path.is_file():
        raise WorkflowError("plan does not exist: %s" % path)
    return run_dir, json.loads(path.read_text(encoding="utf-8"))


def environment_python(project: Path) -> str:
    path = Path(project).resolve() / "config" / "environments.tsv"
    import csv
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row["stage"] in {"orchestrator", "common"} and row["python"].strip():
            return row["python"].strip()
    return sys.executable


def task_command(project: Path, run_id: str, task_id: str) -> list[str]:
    root = Path(project).resolve()
    return [
        environment_python(root), str(root / "Scripts" / "Common" / "run_task.py"),
        "--project", str(root), "--run-id", run_id, "--task", task_id,
    ]


def task_is_implemented(task_id: str, commands: dict) -> bool:
    if task_id == "workflow_summary" or task_id in commands:
        return True
    return any(key.endswith("*") and task_id.startswith(key[:-1]) for key in commands)


def submit(project: Path, run_id: str, dry_run: bool = False) -> dict:
    root = Path(project).resolve()
    run_dir, plan = load_plan(root, run_id)
    current = validate(root)
    if current["status"] != "valid":
        raise WorkflowError("pre-submit preflight failed: %s" % "; ".join(current["errors"]))
    if current["input_signature"] != plan["input_signature"]:
        raise WorkflowError("input signature changed after planning; create a new run-id")
    cfg = load_project(root)
    scheduler = cfg["scheduler"]
    backend = scheduler.get("backend", "local")
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

    job_ids = {}
    submissions = []
    max_parallel = max(1, int(scheduler.get("max_parallel", 2)))
    limited_profiles = set(scheduler.get("limited_profiles") or ["methscan_branch", "dmr", "feature_builder", "trainer"])
    throttle_slots = [None] * max_parallel
    throttle_index = 0
    for item in plan["tasks"]:
        task_id = item["id"]
        snapshot_path = run_dir / "resource_snapshots" / (task_id + ".json")
        snapshot = inspect(root, item["profile"])
        write_json(snapshot_path, snapshot)
        recommendation = snapshot["recommendation"]
        command = task_command(root, run_id, task_id)
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
        }
        if dry_run:
            record["status"] = "dry_run"
            record["job_id"] = "dry_%s" % task_id
            job_ids[task_id] = record["job_id"]
        elif backend == "local":
            env = os.environ.copy()
            env["SCMO_CPUS"] = str(recommendation["cpus"])
            env["SCMO_MEMORY_MB"] = str(recommendation["memory_mb"])
            code = subprocess.call(command, cwd=str(root), env=env)
            if code:
                record["status"] = "failed"
                record["return_code"] = code
                submissions.append(record)
                write_json(run_dir / "submissions.json", submissions)
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
        submissions.append(record)
        write_json(run_dir / "submissions.json", submissions)
    plan["status"] = "dry_run" if dry_run else ("submitted" if backend == "slurm" else "complete")
    plan["job_ids"] = job_ids
    write_json(run_dir / "plan.json", plan)
    return {"status": plan["status"], "run_dir": str(run_dir), "jobs": job_ids}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    result = submit(args.project, args.run_id, args.dry_run)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
