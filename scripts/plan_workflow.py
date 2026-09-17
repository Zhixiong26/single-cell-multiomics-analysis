#!/usr/bin/env python3
"""Build a route-aware task DAG without submitting jobs."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

from _common import WorkflowError, git_commit, load_project, write_json
from validate_project import validate


def task(name: str, profile: str, dependencies: List[str] | None = None, parameters: Dict[str, Any] | None = None) -> dict:
    return {"id": name, "profile": profile, "dependencies": dependencies or [], "parameters": parameters or {}}


def make_tasks(routes: Dict[str, bool], analysis: Dict[str, Any]) -> List[dict]:
    tasks: List[dict] = []
    if routes.get("scanpy"):
        tasks.append(task("scanpy", "scanpy"))
    methscan_tail = None
    if routes.get("methscan_vmr"):
        tasks.extend([
            task("methscan_select_convert", "io_builder"),
            task("methscan_prepare", "serial", ["methscan_select_convert"]),
            task("methscan_filter", "serial", ["methscan_prepare"]),
            task("methscan_smooth", "serial", ["methscan_filter"]),
        ])
        methscan_tail = "methscan_smooth"
        thresholds = analysis.get("methscan", {}).get("vmr_thresholds", [0.01, 0.02, 0.05])
        for threshold in thresholds:
            tasks.append(task("methscan_vmr_%s" % threshold, "methscan_branch", [methscan_tail], {"threshold": threshold}))
    if routes.get("methscan_dmr") and methscan_tail:
        tasks.append(task("methscan_pairwise_dmr", "dmr", [methscan_tail]))
        tasks.append(task("methscan_hypo_heatmaps", "plot", ["methscan_pairwise_dmr"]))
        tasks.append(task("methscan_pooled_dmr", "dmr", [methscan_tail]))
    if routes.get("allcools"):
        deps = ["methscan_filter"] if methscan_tail else []
        tasks.append(task("allcools_features", "feature_builder", deps))
    feature_targets = analysis.get("methylvi", {}).get("feature_targets", [10000, 30000])
    if routes.get("methylvi_allcools"):
        for count in feature_targets:
            tasks.append(task("methylvi_allcools_%s" % count, "trainer", ["allcools_features"], {"features": count}))
    if routes.get("methylvi_vmr") and methscan_tail:
        thresholds = analysis.get("methscan", {}).get("vmr_thresholds", [0.01, 0.02, 0.05])
        for threshold in thresholds:
            vmr_task = "methscan_vmr_%s" % threshold
            for count in feature_targets:
                tasks.append(task("methylvi_vmr_%s_%s" % (threshold, count), "trainer", [vmr_task], {"threshold": threshold, "features": count}))
    if routes.get("methylvi_vmr_dmr") and methscan_tail:
        tasks.append(task("pooled_dmr_prepare", "dmr_prepare", ["methscan_pooled_dmr"]))
        tasks.append(task("pooled_dmr_counts", "feature_builder", ["pooled_dmr_prepare"]))
        previous = "pooled_dmr_counts"
        thresholds = analysis.get("methscan", {}).get("vmr_thresholds", [0.01, 0.02, 0.05])
        for threshold in thresholds:
            for count in feature_targets:
                name = "methylvi_vmr_dmr_%s_%s" % (threshold, count)
                tasks.append(task(name, "trainer", [previous, "methscan_vmr_%s" % threshold], {"threshold": threshold, "features": count}))
                previous = name
    if tasks:
        terminal = {item["id"] for item in tasks}
        for item in tasks:
            terminal.difference_update(item["dependencies"])
        tasks.append(task("workflow_summary", "summary", sorted(terminal)))
    return tasks


def plan(project: Path, routes_option: str, run_id: str) -> dict:
    if not run_id.strip() or "/" in run_id or run_id in {".", ".."}:
        raise WorkflowError("run-id must be a non-empty path-safe name")
    preflight = validate(project)
    if preflight["status"] != "valid":
        raise WorkflowError("preflight failed: %s" % "; ".join(preflight["errors"]))
    available = preflight["routes"]
    if routes_option == "auto":
        selected = dict(available)
    else:
        requested = {item.strip() for item in routes_option.split(",") if item.strip()}
        unknown = requested - set(available)
        blocked = {name for name in requested if not available.get(name)}
        if unknown:
            raise WorkflowError("unknown routes: %s" % ", ".join(sorted(unknown)))
        if blocked:
            raise WorkflowError("routes lack required inputs: %s" % ", ".join(sorted(blocked)))
        selected = {name: name in requested for name in available}
    cfg = load_project(project)
    tasks = make_tasks(selected, cfg["analysis"])
    if not tasks:
        raise WorkflowError("no executable task was selected")
    root = Path(project).resolve()
    run_dir = root / ".workflow" / "runs" / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise WorkflowError("refusing existing non-empty run directory: %s" % run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    template_lock_path = root / ".workflow" / "template-lock.json"
    template_lock = json.loads(template_lock_path.read_text(encoding="utf-8")) if template_lock_path.is_file() else None
    payload = {
        "schema_version": 1, "run_id": run_id, "project": str(root),
        "created_at": datetime.now(timezone.utc).isoformat(), "status": "planned",
        "routes": selected, "input_signature": preflight["input_signature"],
        "git_commit": git_commit(root), "template_lock": template_lock, "tasks": tasks,
    }
    write_json(run_dir / "preflight.json", preflight)
    write_json(run_dir / "plan.json", payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--routes", default="auto")
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    result = plan(args.project, args.routes, args.run_id)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
