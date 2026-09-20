#!/usr/bin/env python3
"""Capture Slurm capacity and recommend a balanced allocation for one task."""

from __future__ import annotations

import argparse
import getpass
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from _common import (WorkflowError, effective_backend, host_role, host_role_source, load_project,
                     parse_memory_mb, write_json)


ACTIVE_STATES = {"idle", "mixed"}


def command_output(command: List[str], fixture: Optional[Path] = None) -> str:
    if fixture:
        return fixture.read_text(encoding="utf-8")
    try:
        return subprocess.check_output(command, stderr=subprocess.STDOUT).decode("utf-8", "replace")
    except (OSError, subprocess.SubprocessError) as exc:
        raise WorkflowError("resource query failed (%s): %s" % (" ".join(command), exc))


def parse_sinfo(text: str) -> List[Dict[str, Any]]:
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("|")
        if len(fields) < 8:
            raise WorkflowError("unexpected sinfo row: %s" % line)
        name, partition, state, total_cpus, cpu_state, real_memory, free_memory, gres = fields[:8]
        allocated, idle, other, total_from_state = [int(value) for value in cpu_state.split("/")]
        rows.append({
            "name": name, "partition": partition.rstrip("*"), "state": state.lower().split("+")[0],
            "cpu_total": int(total_cpus), "cpu_allocated": allocated, "cpu_idle": idle,
            "cpu_other": other, "cpu_state_total": total_from_state,
            "real_memory_mb": int(real_memory), "free_memory_mb": int(free_memory), "gres": gres,
        })
    return rows


def parse_scontrol(text: str) -> Dict[str, Dict[str, Any]]:
    result: Dict[str, Dict[str, Any]] = {}
    for line in text.splitlines():
        values = {}
        for token in line.split():
            if "=" in token:
                key, value = token.split("=", 1)
                values[key] = value
        name = values.get("NodeName")
        if not name:
            continue
        result[name] = {
            "cpu_load": float(values.get("CPULoad", "0") or 0),
            "cpu_allocated": int(values.get("CPUAlloc", "0") or 0),
            "cpu_total": int(values.get("CPUTot", "0") or 0),
            "real_memory_mb": int(values.get("RealMemory", "0") or 0),
            "allocated_memory_mb": int(values.get("AllocMem", "0") or 0),
            "free_memory_mb": int(values.get("FreeMem", "0") or 0),
            "partitions": values.get("Partitions", "").split(","),
            "state_detail": values.get("State", ""),
        }
    return result


def enrich_nodes(sinfo_rows: List[Dict[str, Any]], scontrol_rows: Dict[str, Dict[str, Any]], headroom_mb: int) -> List[Dict[str, Any]]:
    result = []
    for row in sinfo_rows:
        detail = scontrol_rows.get(row["name"], {})
        merged = dict(row)
        merged.update(detail)
        merged["schedulable_cpu"] = max(0, merged["cpu_total"] - merged["cpu_allocated"])
        merged["schedulable_memory_mb"] = max(0, merged["real_memory_mb"] - merged.get("allocated_memory_mb", 0) - headroom_mb)
        merged["observed_free_memory_mb"] = merged.get("free_memory_mb", row["free_memory_mb"])
        result.append(merged)
    return result


def profile_for(scheduler: Dict[str, Any], task: str) -> Dict[str, Any]:
    profiles = scheduler.get("profiles") or {}
    profile = profiles.get(task) or profiles.get("default")
    if not profile:
        raise WorkflowError("scheduler profile is absent for task %s" % task)
    for level in ("floor", "target", "ceiling"):
        if level not in profile:
            raise WorkflowError("resource profile %s lacks %s" % (task, level))
    return profile


def choose(nodes: List[Dict[str, Any]], scheduler: Dict[str, Any], profile: Dict[str, Any]) -> Dict[str, Any]:
    allowed_partitions = set(scheduler.get("partitions") or [])
    allowed_nodes = set(scheduler.get("allow_nodes") or [])
    excluded_nodes = set(scheduler.get("exclude_nodes") or [])
    floor_cpu = int(profile["floor"]["cpus"])
    target_cpu = int(profile["target"]["cpus"])
    ceiling_cpu = int(profile["ceiling"]["cpus"])
    floor_mem = parse_memory_mb(profile["floor"]["memory"])
    target_mem = parse_memory_mb(profile["target"]["memory"])
    ceiling_mem = parse_memory_mb(profile["ceiling"]["memory"])
    candidates = []
    for node in nodes:
        if node["state"] not in ACTIVE_STATES:
            continue
        if allowed_partitions and node["partition"] not in allowed_partitions:
            continue
        if allowed_nodes and node["name"] not in allowed_nodes:
            continue
        if node["name"] in excluded_nodes:
            continue
        if node["schedulable_cpu"] < floor_cpu or node["schedulable_memory_mb"] < floor_mem:
            continue
        cpus = min(target_cpu, ceiling_cpu, node["schedulable_cpu"])
        memory_mb = min(target_mem, ceiling_mem, node["schedulable_memory_mb"])
        if memory_mb < floor_mem:
            continue
        target_fit = int(node["schedulable_cpu"] >= target_cpu and node["schedulable_memory_mb"] >= target_mem)
        waste = (node["schedulable_cpu"] - cpus) + (node["schedulable_memory_mb"] - memory_mb) / 1024.0
        candidates.append((target_fit, -waste, node, cpus, memory_mb))
    if not candidates:
        raise WorkflowError("no active Slurm node satisfies the resource floor")
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    unused_fit, unused_waste, node, cpus, memory_mb = candidates[0]
    return {
        "partition": node["partition"], "reference_node": node["name"],
        "cpus": cpus, "memory_mb": memory_mb, "memory": "%dM" % memory_mb,
        "time": profile["target"].get("time", "1-00:00:00"),
        "pin_node": bool(profile.get("pin_node", False)),
    }


def inspect(project: Path, task: str, sinfo_file: Optional[Path] = None, scontrol_file: Optional[Path] = None, squeue_file: Optional[Path] = None,
            ping_file: Optional[Path] = None, run_id: str = "",
            resolved_backend: Optional[str] = None) -> Dict[str, Any]:
    cfg = load_project(project)
    scheduler = cfg["scheduler"]
    declared = str(scheduler.get("backend") or "local").strip().lower()
    role = host_role(scontrol_file=ping_file) if ping_file is not None else host_role()
    if resolved_backend is None:
        backend, backend_source = effective_backend(scheduler, run_id, role=role)
    else:
        # The caller resolved this through the same rule and holds the
        # acknowledgement in the environment its task will run with, which this
        # process cannot see. It is recorded, not re-derived, so the snapshot and
        # the submission cannot disagree about which backend ran.
        backend = str(resolved_backend).strip().lower()
        backend_source = "declared" if backend == declared else "promoted-from-local"
    snapshot: Dict[str, Any] = {
        "schema_version": 1, "task": task, "backend": backend,
        "backend_declared": declared, "backend_source": backend_source,
        "host_role": role, "host_role_source": host_role_source(),
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    profile = profile_for(scheduler, task)
    if backend == "local":
        limits = scheduler.get("local") or {}
        floor_cpu = int(profile["floor"]["cpus"])
        floor_mem = parse_memory_mb(profile["floor"]["memory"])
        cpus = min(int(profile["target"]["cpus"]), int(limits.get("max_threads", profile["target"]["cpus"])))
        memory_mb = min(parse_memory_mb(profile["target"]["memory"]), parse_memory_mb(limits.get("max_memory", profile["target"]["memory"])))
        if cpus < floor_cpu or memory_mb < floor_mem:
            raise WorkflowError("local limits are below the task resource floor")
        snapshot["recommendation"] = {"cpus": cpus, "memory_mb": memory_mb, "memory": "%dM" % memory_mb}
        snapshot["nodes"] = []
        snapshot["note"] = ("local backend: the task runs on this host, so no cluster is queried and "
                            "the numbers above are the configured limits floored against the profile")
        return snapshot
    if backend != "slurm":
        raise WorkflowError("unsupported scheduler backend: %s" % backend)
    sinfo_text = command_output(["sinfo", "-N", "-h", "-o", "%N|%P|%T|%c|%C|%m|%e|%G"], sinfo_file)
    scontrol_text = command_output(["scontrol", "show", "nodes", "-o"], scontrol_file)
    squeue_command = ["squeue", "-h", "-o", "%i|%u|%T|%P|%R|%C|%m|%M|%l"]
    squeue_text = command_output(squeue_command, squeue_file)
    user_squeue_text = command_output(
        ["squeue", "-h", "-u", getpass.getuser(), "-o", "%i|%u|%T|%P|%R|%C|%m|%M|%l"],
        squeue_file,
    )
    nodes = enrich_nodes(parse_sinfo(sinfo_text), parse_scontrol(scontrol_text), int(scheduler.get("memory_headroom_mb", 4096)))
    snapshot["nodes"] = nodes
    snapshot["queue_global"] = [line for line in squeue_text.splitlines() if line.strip()]
    snapshot["queue_user"] = [line for line in user_squeue_text.splitlines() if line.strip()]
    snapshot["recommendation"] = choose(nodes, scheduler, profile)
    snapshot["note"] = "schedulable memory uses RealMemory-AllocMem-headroom; observed FreeMem is diagnostic only"
    return snapshot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--task", required=True)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--sinfo-file", type=Path)
    parser.add_argument("--scontrol-file", type=Path)
    parser.add_argument("--squeue-file", type=Path)
    args = parser.parse_args()
    result = inspect(args.project, args.task, args.sinfo_file, args.scontrol_file, args.squeue_file)
    if args.json_out:
        write_json(args.json_out, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
