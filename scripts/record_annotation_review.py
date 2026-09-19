#!/usr/bin/env python3
"""Record a reviewed cluster-to-cell-type mapping as the project's annotation profile.

The Scanpy route publishes its cluster labels as proposals: a run leaves
`annotation_status.json` marked pending review, and the cell-type DMR routes refuse those labels.
This tool is the review step in between. It reads the run's own analysis signature and cluster
set, so the profile it writes can only ever describe the run it was made for, and it re-checks
what it is about to write with the same guard the analysis stages use.

Two passes, because the judgement is the reviewer's and this tool has no opinion about it:

  # 1. a worksheet, one row per cluster, seeded with the run's own majority label
  python tools/record_annotation_review.py --project . --run-id run_001 --worksheet review.tsv

  # 2. after editing the cell_type column of that worksheet
  python tools/record_annotation_review.py --project . --run-id run_001 --mapping review.tsv

The profile is written as JSON, which is valid YAML, so recording works whether or not PyYAML is
installed (the loaders accept either). Set `annotation.review_status: approved` in
`config/project.yaml` to let the cell-type DMR routes plan against it.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from _common import WorkflowError, load_project, resolve_path, write_json


DEFAULT_PROFILE = Path("config") / "annotation.yaml"
WORKSHEET_COLUMNS = ("cluster", "proposed_cell_type", "cells", "cell_type", "confidence", "evidence")


def load_json(path: Path, description: str) -> Dict[str, Any]:
    if not path.is_file():
        raise WorkflowError("%s does not exist: %s" % (description, path))
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as exc:
        raise WorkflowError("%s is not valid JSON: %s" % (path, exc))
    if not isinstance(value, dict):
        raise WorkflowError("%s must contain a JSON object: %s" % (description, path))
    return value


def handoff_directory(root: Path, run_id: str) -> Path:
    """Resolve the run's result directory from its own plan rather than by convention.

    `result_dir` is what the planner recorded, so reading it keeps this tool correct even if the
    convention changes; the literal path is only the fallback for a run planned by an older release.
    """
    plan = root / ".workflow" / "runs" / run_id / "plan.json"
    if plan.is_file():
        try:
            recorded = json.loads(plan.read_text(encoding="utf-8")).get("result_dir")
        except ValueError:
            recorded = None
        if recorded:
            return Path(str(recorded)).resolve()
    return (root / "Results" / "runs" / run_id).resolve()


def read_annotation_table(path: Path, cluster_column: str) -> List[Dict[str, str]]:
    if not path.is_file():
        raise WorkflowError(
            "the scanpy task published no annotation table at %s; run the scanpy task for this "
            "run before reviewing its clusters" % path
        )
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = list(reader.fieldnames or ())
        missing = [name for name in ("cell_id", "cell_type", cluster_column) if name not in columns]
        if missing:
            raise WorkflowError(
                "%s lacks column(s) %s; the published handoff has %s. Pass --cluster-column to name "
                "the cluster column explicitly."
                % (path, ", ".join(missing), ", ".join(columns) or "none")
            )
        rows = list(reader)
    if not rows:
        raise WorkflowError("the published annotation table has no cells: %s" % path)
    return rows


def summarise_clusters(rows: List[Dict[str, str]], cluster_column: str) -> Dict[str, Dict[str, Any]]:
    """Count cells and tally labels per cluster, with the majority label as the proposal."""
    clusters: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        cluster = (row.get(cluster_column) or "").strip()
        if not cluster:
            raise WorkflowError("the annotation table has an empty %s value" % cluster_column)
        label = (row.get("cell_type") or "").strip()
        entry = clusters.setdefault(cluster, {"cells": 0, "labels": {}})
        entry["cells"] += 1
        entry["labels"][label] = entry["labels"].get(label, 0) + 1
    for entry in clusters.values():
        # Sorted by descending count then label so an exact tie is resolved the same way every run.
        entry["proposed"] = sorted(entry["labels"].items(), key=lambda item: (-item[1], item[0]))[0][0]
    return clusters


def read_worksheet(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        raise WorkflowError("no worksheet to record: %s does not exist" % path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "cluster" not in (reader.fieldnames or []) or "cell_type" not in (reader.fieldnames or []):
            raise WorkflowError("%s needs cluster and cell_type columns" % path)
        return list(reader)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", help="run whose published annotation handoff is under review")
    parser.add_argument("--annotation-status", type=Path,
                        help="point at a handoff status file directly instead of --run-id")
    parser.add_argument("--cluster-column", default="leiden",
                        help="annotation table column holding the clusters (default: leiden)")
    parser.add_argument("--worksheet", type=Path, help="write a reviewer worksheet and stop")
    parser.add_argument("--mapping", type=Path, help="read a completed worksheet and record it")
    args = parser.parse_args()
    if bool(args.worksheet) == bool(args.mapping):
        raise WorkflowError("pass exactly one of --worksheet (prepare the review) or "
                            "--mapping (record it)")
    if not args.run_id and not args.annotation_status:
        raise WorkflowError("pass --run-id, or --annotation-status to name a handoff directly")

    root = args.project.resolve()
    cfg = load_project(root)
    annotation = cfg.get("annotation") or {}
    if args.annotation_status:
        status_path = args.annotation_status.resolve()
        table_path = status_path.with_name("annotation.tsv")
    else:
        directory = handoff_directory(root, str(args.run_id))
        status_path = directory / "annotation_status.json"
        table_path = directory / "annotation.tsv"
    status = load_json(status_path, "the scanpy annotation handoff")
    rows = read_annotation_table(table_path, args.cluster_column)
    clusters = summarise_clusters(rows, args.cluster_column)

    signature = status.get("analysis_signature")
    if not signature:
        raise WorkflowError(
            "the handoff at %s records no analysis_signature, so a profile recorded from it could "
            "not be checked against the run it describes" % status_path
        )
    run_kind = status.get("run_kind")

    if args.worksheet:
        destination = args.worksheet.resolve()
        if destination.exists():
            # Overwriting would destroy a review in progress; a worksheet is the reviewer's work.
            raise WorkflowError("%s already exists; move it aside or choose another path" % destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=WORKSHEET_COLUMNS, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            for cluster in sorted(clusters):
                entry = clusters[cluster]
                writer.writerow({
                    "cluster": cluster, "proposed_cell_type": entry["proposed"], "cells": entry["cells"],
                    "cell_type": entry["proposed"], "confidence": "reviewed", "evidence": "",
                })
        minor = {cluster: entry["labels"] for cluster, entry in clusters.items()
                 if len(entry["labels"]) > 1}
        print("worksheet: %s (%d cluster(s), %d cell(s) from the %s run)"
              % (destination, len(clusters), len(rows), run_kind or "scanpy"))
        print("cell_type is pre-filled with each cluster's majority label; correct every row you "
              "disagree with, then record the review with --mapping.", flush=True)
        if minor:
            # A cluster whose cells disagree is the case most likely to be mislabelled by majority vote.
            print("clusters whose cells carry more than one label (check these first):", file=sys.stderr)
            for cluster in sorted(minor):
                tally = ", ".join("%s=%d" % item for item in
                                  sorted(minor[cluster].items(), key=lambda item: (-item[1], item[0])))
                print("  %s: %s" % (cluster, tally), file=sys.stderr)
        return 0

    # Local import on purpose, and used only in this function: the path to the guard's module
    # depends on the project root, which is not known until argparse has run. Module-level helpers
    # above must never reference it.
    sys.path.insert(0, str(root / "Scripts" / "Common"))
    from annotation import PLACEHOLDER_LABELS, evaluate_annotation_profile  # noqa: E402

    mapping: Dict[str, str] = {}
    confidence: Dict[str, str] = {}
    evidence: Dict[str, str] = {}
    for row in read_worksheet(args.mapping.resolve()):
        cluster = (row.get("cluster") or "").strip()
        if not cluster:
            raise WorkflowError("the worksheet has a row with no cluster")
        if cluster in mapping:
            raise WorkflowError("the worksheet lists cluster %s twice" % cluster)
        label = (row.get("cell_type") or "").strip()
        if not label:
            raise WorkflowError("cluster %s has no cell_type in the worksheet" % cluster)
        if label.lower() in PLACEHOLDER_LABELS:
            # Recording one would hand the cell-type DMR stages a "cell type" that they skip or pool
            # without saying so (Scripts/Methscan/07_methdiff_celltype.py), so refuse it here where
            # the reviewer can still fix it.
            raise WorkflowError(
                "cluster %s is mapped to the placeholder %r, which is not a cell type; give it a "
                "real cell type or leave the profile unrecorded" % (cluster, label)
            )
        mapping[cluster] = label
        confidence[cluster] = (row.get("confidence") or "").strip() or "reviewed"
        evidence[cluster] = (row.get("evidence") or "").strip()

    observed = sorted(clusters)
    recorded = sorted(mapping)
    if recorded != observed:
        absent = [cluster for cluster in observed if cluster not in mapping]
        unknown = [cluster for cluster in recorded if cluster not in clusters]
        raise WorkflowError(
            "the worksheet does not describe this run: missing %s, not in the run %s. Re-run with "
            "--worksheet to review the current cluster set."
            % (absent or "nothing", unknown or "nothing")
        )

    profile: Dict[str, Any] = {
        "analysis_signature": signature,
        "expected_clusters": observed,
        "annotations": {
            cluster: {"cell_type": mapping[cluster], "confidence": confidence[cluster],
                      "evidence": evidence[cluster]}
            for cluster in observed
        },
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "recorded_from": str(status_path),
    }
    _, verdict, _template = evaluate_annotation_profile(observed, signature, profile)
    if verdict["status"] != "reviewed":
        # Unreachable while the checks above hold, and kept because this is an approval gate: the
        # profile must pass the same guard a baseline run will hold it to before it reaches disk.
        raise WorkflowError("refusing to record a profile the annotation guard rejects (%s)"
                            % verdict.get("reason"))

    configured = resolve_path(root, annotation.get("profile"))
    destination = configured or (root / DEFAULT_PROFILE)
    write_json(destination, profile)
    print("recorded %d cluster(s) with analysis signature %s\n  %s"
          % (len(observed), signature, destination))

    if configured is None:
        print("\nconfig/project.yaml does not set annotation.profile, so the runners will not look "
              "at that file. Add the key to its existing annotation block:\n"
              "  annotation:\n    profile: config/annotation.yaml", file=sys.stderr)
    if str(annotation.get("review_status") or "") != "approved":
        print("\nannotation.review_status in config/project.yaml is %r; the cell-type DMR routes "
              "plan only when it is 'approved'. Set it once this mapping has been checked."
              % annotation.get("review_status"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
