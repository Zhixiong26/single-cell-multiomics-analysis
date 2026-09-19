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
    WorkflowError, load_environments, load_project, load_samples, parse_memory_mb,
    plan_data_sources, project_files, resolve_path, sha256_file, signature,
    stored_data_sources, write_json,
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


def validate(project: Path, require_paths: bool = True, mode: str = "quick", workers: int = 1,
             required_stages: set[str] | None = None,
             selected_routes: set[str] | None = None) -> dict:
    if mode not in {"quick", "full"}:
        raise WorkflowError("validation mode must be quick or full")
    files = project_files(project)
    cfg = load_project(project)
    samples = load_samples(project, require_paths=require_paths)
    active = [row for row in samples if row["include"]]
    errors, warnings = [], []
    if int(cfg.get("schema_version", 1)) == 1:
        warnings.append("schema v1 is deprecated; migrate the project configuration to schema v2")
    all_rna_samples = [row for row in active if row["rna_path"]]
    all_allc_samples = [row for row in active if row["allc_root"]]
    use_rna = selected_routes is None or "scanpy" in selected_routes
    use_allc = selected_routes is None or bool(selected_routes - {"scanpy"})
    rna_samples = all_rna_samples if use_rna else []
    allc_samples = all_allc_samples if use_allc else []
    allc_inventory, allc_cell_ids, validation_jobs = [], [], []
    context = str(cfg.get("analysis", {}).get("allcools", {}).get("mc_context", "CGN"))
    thresholds = cfg.get("analysis", {}).get("methscan", {}).get("vmr_thresholds", [0.01, 0.02, 0.05])
    needs_thresholds = selected_routes is None or bool(selected_routes & {
        "methscan_vmr", "methscan_dmr", "methylvi_vmr", "methylvi_vmr_dmr",
    })
    if needs_thresholds and (not isinstance(thresholds, list) or not thresholds):
        errors.append("analysis.methscan.vmr_thresholds must be a non-empty list")
    elif needs_thresholds:
        try:
            numeric_thresholds = [float(value) for value in thresholds]
            if any(value <= 0 or value > 1 for value in numeric_thresholds):
                errors.append("analysis.methscan.vmr_thresholds values must satisfy 0 < value <= 1")
            if len(set(numeric_thresholds)) != len(numeric_thresholds):
                errors.append("analysis.methscan.vmr_thresholds values must be unique")
        except (TypeError, ValueError):
            errors.append("analysis.methscan.vmr_thresholds values must be numeric")
    feature_targets = cfg.get("analysis", {}).get("methylvi", {}).get("feature_targets", [10000, 30000])
    needs_features = selected_routes is None or bool(selected_routes & {
        "allcools", "methylvi_allcools", "methylvi_vmr", "methylvi_vmr_dmr",
    })
    if needs_features and (not isinstance(feature_targets, list) or not feature_targets):
        errors.append("analysis.methylvi.feature_targets must be a non-empty list")
    elif needs_features:
        valid_targets = all(isinstance(value, int) and not isinstance(value, bool) and value > 0
                            for value in feature_targets)
        if not valid_targets:
            errors.append("analysis.methylvi.feature_targets values must be positive integers")
        elif len(set(feature_targets)) != len(feature_targets):
            errors.append("analysis.methylvi.feature_targets values must be unique")
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

    # Declared inputs are checked where they will be read, not where they were
    # declared: a link that resolves on the login host and dangles on the compute
    # host is the failure this catches, and it names the link rather than the stage.
    data_plans = plan_data_sources(files["root"], stored_data_sources(cfg))
    data_validation = []
    for plan in data_plans:
        entry = {"name": plan["name"], "destination": plan["destination"],
                 "state": plan["state"], "origin": plan["url"] or plan["target"]}
        if plan["state"] == "conflict":
            errors.append("Data/%s exists but is not the entry data.sources declares (%s)"
                          % (plan["name"], entry["origin"]))
        elif plan["state"] == "dangling" and require_paths:
            errors.append("Data/%s does not resolve: %s is not reachable from here, so every stage would "
                          "fail on it. Run `python tools/link_data.py --project . --execute` where both are "
                          "visible, or validate from the execution context" % (plan["name"], entry["origin"]))
        elif plan["state"] == "download" and require_paths:
            errors.append("Data/%s is declared as a download and has not been fetched; run "
                          "`python tools/link_data.py --project . --execute`" % plan["name"])
        elif plan["state"] == "present" and mode == "full":
            # Deferred to full mode: the digest is the only thing that can tell a
            # finished transfer from a truncated one, and it reads the whole file.
            observed = sha256_file(Path(plan["destination"]))
            if observed.lower() != plan["sha256"]:
                errors.append("Data/%s has sha256 %s, not the declared %s"
                              % (plan["name"], observed, plan["sha256"]))
            entry["sha256"] = observed
        data_validation.append(entry)

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
    if required_stages is None:
        required_stages = {"orchestrator", "scanpy_allcools"}
        if allc_samples:
            required_stages.update({"methscan", "methylvi"})
    for env in environment_rows:
        declared_required = env["required"].strip().lower() in {"1", "true", "yes"}
        required = env["stage"].strip() in required_stages
        executable = resolve_path(files["root"], env["executable"])
        python = resolve_path(files["root"], env["python"])
        item = {"stage": env["stage"], "required": required, "declared_required": declared_required,
                "python": str(python) if python else "",
                "executable": str(executable) if executable else ""}
        if required and require_paths and (python is None or not python.is_file()):
            errors.append("required Python is absent for %s: %s" % (env["stage"], python))
        if required and require_paths and (executable is None or not executable.is_file()):
            errors.append("required executable is absent for %s: %s" % (env["stage"], executable))
            item["status"] = "missing"
        else:
            item["status"] = "present" if executable and executable.exists() else "optional_absent"
        if required and env["version_command"].strip() and item["status"] == "present":
            try:
                output = subprocess.check_output(shlex.split(env["version_command"]), stderr=subprocess.STDOUT, timeout=30)
                item["version"] = output.decode("utf-8", "replace").strip().splitlines()[0]
            except (OSError, subprocess.SubprocessError) as exc:
                errors.append("version check failed for %s: %s" % (env["stage"], exc))
        env_results.append(item)
    declared_stages = {row["stage"].strip() for row in environment_rows}
    missing_stages = sorted(required_stages - declared_stages)
    if missing_stages:
        errors.append("required environment stages are absent: %s; run tools/bootstrap_environments.py --project PROJECT --execute" % ", ".join(missing_stages))

    scheduler = cfg.get("scheduler") or {}
    backend = scheduler.get("backend", "local")
    if backend not in {"local", "slurm"}:
        errors.append("scheduler.backend must be local or slurm")
    if backend == "slurm" and not scheduler.get("partitions"):
        errors.append("Slurm backend requires a non-empty scheduler.partitions allow-list")
    profiles = scheduler.get("profiles") or {}
    required_profiles = {"scanpy", "io_builder", "serial", "methscan_branch", "dmr",
                         "dmr_prepare", "feature_builder", "trainer", "plot", "summary"}
    for name in sorted(required_profiles):
        profile = profiles.get(name) or profiles.get("default")
        if not profile or any(level not in profile for level in ("floor", "target", "ceiling")):
            errors.append("scheduler resource profile %s lacks floor/target/ceiling" % name)
    for name, profile in sorted(profiles.items()):
        # A profile that is not monotonic cannot mean what its author intended: inspect_resources
        # clamps the target down to the ceiling and refuses any node below the floor, so a floor
        # above the target asks for more than the job will ever request. The `dmr` floor once sat
        # below the parallelism its own stage declares -- the same class of mistake -- so the
        # ordering is checked here rather than discovered as a failed task on a busy node.
        if not isinstance(profile, dict) or any(level not in profile for level in ("floor", "target", "ceiling")):
            continue  # already reported above when the profile is one of the required ones
        levels = ("floor", "target", "ceiling")
        try:
            cpus = {level: int(profile[level]["cpus"]) for level in levels}
            memory = {level: parse_memory_mb(profile[level]["memory"]) for level in levels}
        except (KeyError, TypeError, ValueError, WorkflowError) as exc:
            errors.append("scheduler resource profile %s has an unreadable cpus/memory value: %s"
                          % (name, exc))
            continue
        for field, values in (("cpus", cpus), ("memory", memory)):
            if not values["floor"] <= values["target"] <= values["ceiling"]:
                errors.append(
                    "scheduler resource profile %s must satisfy floor <= target <= ceiling for %s; "
                    "it declares %s"
                    % (name, field, " <= ".join(str(values[level]) for level in levels)))

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
        "scanpy": bool(all_rna_samples), "methscan_vmr": bool(all_allc_samples), "methscan_dmr": dmr_ready,
        "allcools": bool(all_allc_samples), "methylvi_allcools": bool(all_allc_samples),
        "methylvi_vmr": bool(all_allc_samples), "methylvi_vmr_dmr": dmr_ready,
    }
    if allc_samples and not annotation_available:
        warnings.append("cell-type DMR routes are disabled until an annotation table is provided")

    rna_inventory = []
    for row in rna_samples:
        path = Path(row["rna_path"])
        stat = path.stat()
        rna_inventory.append({"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    code_hashes, code_suffixes = {}, {".py", ".sh", ".sbatch"}
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
        "selected_routes": sorted(selected_routes) if selected_routes is not None else None,
        "counts": {"samples": len(active), "rna_samples": len(rna_samples),
                   "allc_samples": len(allc_samples), "allc_cells": len(allc_inventory)},
        "routes": routes, "references": reference_hashes, "environments": env_results,
        "allc_validation": sorted(allc_validation, key=lambda item: item["path"]),
        # Recorded so a reader can see where each input actually came from. Not
        # signed: the declared sources are already part of the input signature, and
        # a link's transient state would otherwise churn that signature for nothing.
        "data": data_validation,
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
