#!/usr/bin/env python3
"""Read-only quick/full preflight for inputs, references, environments, and routes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from _common import (
    WorkflowError, load_environments, load_project, load_samples, project_files,
    resolve_path, sha256_file, signature, write_json,
)


def context_matches(observed: str, expected: str) -> bool:
    observed, expected = observed.upper(), expected.upper()
    return observed.startswith(expected[:-1]) if expected.endswith("N") else observed == expected


def validate_allc_record(path: Path, context: str = "CGN", full: bool = False) -> dict:
    opener = gzip.open if path.name.endswith(".gz") else open
    saw_context, records, current_chrom, last_position, seen_chroms = False, 0, None, None, set()
    with opener(str(path), "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise WorkflowError("%s:%d has fewer than 6 ALLC columns" % (path, line_number))
            try:
                position, mc, cov = int(fields[1]), int(fields[4]), int(fields[5])
            except ValueError:
                if line_number == 1:
                    continue
                raise WorkflowError("%s:%d has invalid integer fields" % (path, line_number))
            if position < 1 or mc < 0 or cov < 0 or mc > cov:
                raise WorkflowError("%s:%d has invalid ALLC counts" % (path, line_number))
            if full:
                if fields[0] != current_chrom:
                    if fields[0] in seen_chroms:
                        raise WorkflowError("%s:%d repeats a completed chromosome block" % (path, line_number))
                    seen_chroms.add(fields[0])
                    current_chrom, last_position = fields[0], None
                if last_position is not None and position < last_position:
                    raise WorkflowError("%s:%d is not coordinate sorted" % (path, line_number))
                last_position = position
            saw_context = saw_context or context_matches(fields[3], context)
            records += 1
            if not full and records >= 100 and saw_context:
                break
    if not saw_context:
        raise WorkflowError("ALLC preflight found no %s context in checked records: %s" % (context, path))
    return {"path": str(path.resolve()), "records_checked": records,
            "mode": "full" if full else "quick"}


def allc_cell_id(path: Path, row: dict) -> str:
    regex = row.get("allc_cell_id_regex", "").strip()
    replacement = row.get("allc_cell_id_replacement", "")
    if regex:
        match = re.search(regex, path.name)
        if not match:
            raise WorkflowError("ALLC filename does not match configured cell-ID regex: %s" % path)
        value = match.expand(replacement) if replacement else (match.groupdict().get("cell_id") or match.group(0))
    else:
        value = path.name
        for suffix in (".allc.tsv.gz", "_allc.gz", ".allc.gz"):
            if value.endswith(suffix):
                value = value[:-len(suffix)]
                break
    if not value:
        raise WorkflowError("ALLC filename produced an empty cell ID: %s" % path)
    prefix = row["cell_id_prefix"]
    return value if not prefix or value.startswith(prefix + "_") else prefix + "_" + value


def validate(project: Path, require_paths: bool = True, mode: str = "quick", workers: int = 1) -> dict:
    if mode not in {"quick", "full"}:
        raise WorkflowError("validation mode must be quick or full")
    files = project_files(project)
    cfg = load_project(project)
    samples = load_samples(project, require_paths=require_paths)
    active = [row for row in samples if row["include"]]
    errors, warnings = [], []
    rna_samples = [row for row in active if row["rna_path"]]
    allc_samples = [row for row in active if row["allc_root"]]
    allc_inventory, allc_cell_ids, validation_jobs = [], [], []
    context = str(cfg.get("analysis", {}).get("allcools", {}).get("mc_context", "CGN"))
    for row in allc_samples:
        root = Path(row["allc_root"])
        pattern = row["allc_glob"].strip() or "**/*.allc.tsv.gz"
        paths = sorted(path for path in root.glob(pattern) if path.is_file()) if root.exists() else []
        if require_paths and not paths:
            errors.append("sample %s has no ALLC matching %s" % (row["sample_id"], pattern))
        for path in paths:
            if not Path(str(path) + ".tbi").is_file():
                errors.append("indexed ALLC is missing .tbi: %s" % path)
            try:
                allc_cell_ids.append(allc_cell_id(path, row))
            except WorkflowError as exc:
                errors.append(str(exc))
            stat = path.stat()
            allc_inventory.append({"path": str(path.resolve()), "size": stat.st_size,
                                   "mtime_ns": stat.st_mtime_ns, "sample_id": row["sample_id"]})
            if require_paths:
                validation_jobs.append(path)
    allc_validation = []
    if validation_jobs:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = {pool.submit(validate_allc_record, path, context, mode == "full"): path
                       for path in validation_jobs}
            for future in as_completed(futures):
                try:
                    allc_validation.append(future.result())
                except (OSError, EOFError, WorkflowError) as exc:
                    errors.append(str(exc))
    if len(set(allc_cell_ids)) != len(allc_cell_ids):
        errors.append("ALLC-derived cell IDs must be unique across included samples")

    references = cfg.get("references") or {}
    reference_hashes = {}
    for key in ("chrom_sizes", "blacklist", "tss_bed"):
        value = references.get(key)
        if not value:
            continue
        path = resolve_path(files["root"], value)
        if require_paths and (path is None or not path.is_file()):
            errors.append("reference %s does not exist: %s" % (key, path))
        elif path and path.is_file():
            reference_hashes[key] = sha256_file(path)
            expected = references.get(key + "_sha256")
            if expected and expected.lower() != reference_hashes[key].lower():
                errors.append("reference checksum mismatch for %s" % key)
    if allc_samples:
        for key in ("chrom_sizes", "blacklist"):
            if not references.get(key):
                errors.append("methylation routes require reference %s" % key)

    annotation = cfg.get("annotation") or {}
    legacy = annotation.get("path")
    if legacy and not annotation.get("table"):
        warnings.append("annotation.path is deprecated; migrate it to annotation.table")
    annotation_path = resolve_path(files["root"], annotation.get("table") or legacy)
    annotation_profile = resolve_path(files["root"], annotation.get("profile"))
    if annotation_path and require_paths and not annotation_path.is_file():
        errors.append("annotation does not exist: %s" % annotation_path)
    if annotation_profile and require_paths and not annotation_profile.is_file():
        errors.append("annotation profile does not exist: %s" % annotation_profile)
    annotation_available = bool(annotation_path and (annotation_path.is_file() or not require_paths))
    annotation_rows, annotation_hash = [], None
    cell_column = annotation.get("cell_id_column", "cell_id")
    type_column = annotation.get("cell_type_column", "cell_type")
    if annotation_path and annotation_path.is_file():
        annotation_hash = sha256_file(annotation_path)
        with annotation_path.open("r", encoding="utf-8", newline="") as handle:
            annotation_rows = list(csv.DictReader(handle, delimiter="\t"))
        if not annotation_rows or cell_column not in annotation_rows[0] or type_column not in annotation_rows[0]:
            errors.append("annotation must contain %s and %s" % (cell_column, type_column))
        elif len({row[cell_column] for row in annotation_rows}) != len(annotation_rows):
            errors.append("annotation cell IDs must be unique")
        elif any(not row[cell_column].strip() or not row[type_column].strip() for row in annotation_rows):
            errors.append("annotation cell IDs and cell types must be non-empty")

    env_results = []
    environment_rows = load_environments(project)
    for env in environment_rows:
        required = env["required"].strip().lower() in {"1", "true", "yes"}
        executable = resolve_path(files["root"], env["executable"])
        python = resolve_path(files["root"], env["python"])
        item = {"stage": env["stage"], "required": required, "python": str(python) if python else "",
                "executable": str(executable) if executable else ""}
        if required and require_paths and (python is None or not python.is_file()):
            errors.append("required Python is absent for %s: %s" % (env["stage"], python))
        if required and require_paths and (executable is None or not executable.is_file()):
            errors.append("required executable is absent for %s: %s" % (env["stage"], executable))
            item["status"] = "missing"
        else:
            item["status"] = "present" if executable and executable.exists() else "optional_absent"
        if env["version_command"].strip() and item["status"] == "present":
            try:
                output = subprocess.check_output(shlex.split(env["version_command"]), stderr=subprocess.STDOUT, timeout=30)
                item["version"] = output.decode("utf-8", "replace").strip().splitlines()[0]
            except (OSError, subprocess.SubprocessError) as exc:
                errors.append("version check failed for %s: %s" % (env["stage"], exc))
        env_results.append(item)
    declared_stages = {row["stage"].strip() for row in environment_rows}
    required_stages = {"orchestrator", "scanpy_allcools"}
    if allc_samples:
        required_stages.update({"methscan", "methylvi"})
    missing_stages = sorted(required_stages - declared_stages)
    if missing_stages:
        errors.append("required environment stages are absent: %s; run tools/bootstrap_environments.py --project PROJECT --execute" % ", ".join(missing_stages))

    minimum = int(cfg.get("analysis", {}).get("methscan", {}).get("min_cells", 6))
    selected_annotation = annotation_rows
    if allc_cell_ids and annotation_rows and cell_column in annotation_rows[0] and type_column in annotation_rows[0]:
        allc_set = set(allc_cell_ids)
        selected_annotation = [row for row in annotation_rows if row[cell_column] in allc_set]
        if not selected_annotation:
            warnings.append("annotation contains no cell IDs matching included ALLC inputs")
    annotated_types = {}
    for row in selected_annotation:
        label = row.get(type_column, "").strip()
        if label:
            annotated_types[label] = annotated_types.get(label, 0) + 1
    eligible_types = {label: count for label, count in annotated_types.items() if count >= minimum}
    forbidden = sorted(set(annotated_types) & {"NA", "Unassigned", "requires_review"})
    approved = annotation.get("review_status", "unreviewed") == "approved"
    dmr_ready = bool(allc_samples and annotation_available and approved and len(eligible_types) >= 2 and not forbidden)
    if allc_samples and annotation_available and not dmr_ready:
        warnings.append("cell-type DMR routes require approved annotation, no placeholder labels, and at least two cell types with >= %d matching cells" % minimum)
    routes = {
        "scanpy": bool(rna_samples), "methscan_vmr": bool(allc_samples), "methscan_dmr": dmr_ready,
        "allcools": bool(allc_samples), "methylvi_allcools": bool(allc_samples),
        "methylvi_vmr": bool(allc_samples), "methylvi_vmr_dmr": dmr_ready,
    }
    if allc_samples and not annotation_available:
        warnings.append("cell-type DMR routes are disabled until an annotation table is provided")

    rna_inventory = []
    for row in rna_samples:
        path = Path(row["rna_path"])
        stat = path.stat()
        rna_inventory.append({"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    code_hashes, code_suffixes = {}, {".py", ".sh", ".sbatch", ".ipynb"}
    for code_root in (files["root"] / "Scripts", files["root"] / "tools"):
        if code_root.exists():
            for path in sorted(item for item in code_root.rglob("*") if item.is_file() and item.suffix in code_suffixes):
                code_hashes[str(path.relative_to(files["root"]))] = sha256_file(path)
    spec_root = files["root"] / "environment-specs"
    if spec_root.exists():
        for path in sorted(spec_root.glob("*.yaml")):
            code_hashes[str(path.relative_to(files["root"]))] = sha256_file(path)
    workflow_entry = files["root"] / "workflow.py"
    if workflow_entry.is_file():
        code_hashes[str(workflow_entry.relative_to(files["root"]))] = sha256_file(workflow_entry)
    code_signature = signature(code_hashes)
    payload = {
        "status": "invalid" if errors else "valid", "schema_version": 2,
        "project": str(files["root"]), "validation_mode": mode,
        "counts": {"samples": len(active), "rna_samples": len(rna_samples),
                   "allc_samples": len(allc_samples), "allc_cells": len(allc_inventory)},
        "routes": routes, "references": reference_hashes, "environments": env_results,
        "allc_validation": sorted(allc_validation, key=lambda item: item["path"]),
        "code_signature": code_signature, "errors": errors, "warnings": warnings,
    }
    payload["input_signature"] = signature({
        "project": {key: value for key, value in cfg.items() if key not in {"project_root", "config_files"}},
        "samples": active, "rna_inventory": rna_inventory, "allc_inventory": allc_inventory,
        "allc_cell_ids": allc_cell_ids, "annotation_sha256": annotation_hash,
        "code_signature": code_signature, "references": reference_hashes,
        "environments": env_results, "routes": routes,
    })
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--allow-missing-paths", action="store_true")
    parser.add_argument("--mode", choices=("quick", "full"), default="quick")
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    args = parser.parse_args()
    result = validate(args.project, require_paths=not args.allow_missing_paths,
                      mode=args.mode, workers=args.workers)
    if args.json_out:
        write_json(args.json_out, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "valid" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
