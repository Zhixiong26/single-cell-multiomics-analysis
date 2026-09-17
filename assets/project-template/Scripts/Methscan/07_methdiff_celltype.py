#!/usr/bin/env python3
"""Run pairwise MethSCAn differential methylation between RNA cell types."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


PRIMARY_CHROM_RE = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methscan", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Filtered and smoothed MethSCAn data directory")
    parser.add_argument("--cell-metadata", type=Path, required=True,
                        help="Manifest containing cell_id and rna_cell_type")
    parser.add_argument("--expected-cell-ids", type=Path, required=True,
                        help="Filtered column_header.txt")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-cells", type=int, default=6)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=2,
                        help="Concurrent pairwise comparisons")
    parser.add_argument("--bandwidth", type=int, default=2000)
    parser.add_argument("--stepsize", type=int, default=1000)
    parser.add_argument("--threshold", type=float, default=0.02)
    parser.add_argument("--bridge-gaps", type=int, default=0)
    parser.add_argument("--resume", action="store_true",
                        help="Reuse validated complete comparisons in a partial output directory")
    parser.add_argument(
        "--pool-samples-as",
        default=None,
        help="Pool cells across sample_id values and run one cell-type analysis under this label",
    )
    return parser.parse_args()


def sanitize(label):
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(label)).strip("_")
    if not value:
        raise ValueError("Cannot sanitize empty cell type")
    return value


def read_metadata(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"cell_id", "sample_id", "rna_cell_type"}.issubset(rows[0]):
        raise ValueError("Cell metadata must contain cell_id, sample_id and rna_cell_type")
    cells = {}
    for row in rows:
        cell_id = row["cell_id"]
        sample = (row.get("sample_id") or "").strip()
        cell_type = (row.get("rna_cell_type") or "").strip()
        if not cell_id or not sample or not cell_type or cell_type == "NA":
            continue
        value = (sample, cell_type)
        if cell_id in cells and cells[cell_id] != value:
            raise ValueError("Conflicting sample/cell type for %s" % cell_id)
        cells[cell_id] = value
    return cells


def read_filtered_ids(path):
    with path.open() as handle:
        ids = [line.strip() for line in handle if line.strip()]
    if len(ids) != len(set(ids)):
        raise ValueError("Filtered cell IDs are not unique")
    return ids


def build_primary_data_view(source_dir, output_dir):
    """Create a no-copy hard-link view containing chr1-22/X/Y only."""
    view = output_dir / "methscan_input_primary"
    if view.is_dir():
        manifest = view / "primary_view.json"
        if not manifest.is_file():
            raise FileNotFoundError(manifest)
        recorded = json.loads(manifest.read_text())
        if Path(recorded.get("source_data_dir", "")).resolve() != source_dir.resolve():
            raise ValueError("Existing primary view points to a different source data directory")
        expected = {"chr%d.npz" % value for value in range(1, 23)} | {"chrX.npz", "chrY.npz"}
        found = {path.name for path in view.glob("*.npz")}
        smooth_found = {path.stem + ".npz" for path in (view / "smoothed").glob("*.csv")}
        if found != expected or smooth_found != expected or not (view / "column_header.txt").is_file():
            raise ValueError("Existing primary chromosome view is incomplete")
        return view
    smooth_view = view / "smoothed"
    smooth_view.mkdir(parents=True)
    linked = []
    for name in ("column_header.txt", "cell_stats.csv", "run_info.txt"):
        source = source_dir / name
        if source.is_file():
            os.link(source, view / name)
            linked.append(name)
    if not (view / "column_header.txt").is_file():
        raise FileNotFoundError(source_dir / "column_header.txt")
    for source in sorted(source_dir.glob("*.npz")):
        chrom = source.stem
        if not PRIMARY_CHROM_RE.fullmatch(chrom):
            continue
        smooth_source = source_dir / "smoothed" / (chrom + ".csv")
        if not smooth_source.is_file():
            raise FileNotFoundError(smooth_source)
        os.link(source, view / source.name)
        os.link(smooth_source, smooth_view / smooth_source.name)
        linked.append(source.name)
    expected = {"chr%d.npz" % value for value in range(1, 23)} | {"chrX.npz", "chrY.npz"}
    found = {name for name in linked if name.endswith(".npz")}
    if found != expected:
        raise ValueError("Primary chromosome set mismatch: missing=%s extra=%s" %
                         (sorted(expected - found), sorted(found - expected)))
    (view / "primary_view.json").write_text(json.dumps({
        "source_data_dir": str(source_dir.resolve()),
        "link_type": "hard_link_no_copy",
        "chromosomes": sorted(path.stem for path in view.glob("*.npz")),
    }, indent=2, sort_keys=True) + "\n")
    return view


def validate_dmr(path):
    if not path.is_file():
        return False
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                return False
            if not fields[1].isdigit() or not fields[2].isdigit() or int(fields[2]) <= int(fields[1]):
                return False
    return True


def is_known_fdr_zero_failure(path):
    if not path.is_file():
        return False
    text = path.read_text(errors="replace")
    return (
        'calc_fdr(output_final[11] == "real")' in text
        and "ZeroDivisionError: division by zero" in text
    )


def run_comparison(task):
    args, comparison = task
    sample, label_a, label_b, cells_a, cells_b = comparison
    pair = "%s_vs_%s" % (sanitize(label_a), sanitize(label_b))
    pair_dir = args.output_dir / "samples" / sanitize(sample) / "comparisons" / pair
    group_file = pair_dir / "cell_groups.csv"
    output_file = pair_dir / "DMRs.bed"
    log_file = pair_dir / "methscan_diff.log"
    status_file = pair_dir / "status.json"
    pair_dir.mkdir(parents=True, exist_ok=True)

    if args.resume and status_file.is_file():
        try:
            previous = json.loads(status_file.read_text())
        except (OSError, json.JSONDecodeError):
            previous = {}
        identity_matches = all([
            previous.get("sample_id") == sample,
            previous.get("comparison") == pair,
            previous.get("cell_type_a") == label_a,
            previous.get("cell_type_b") == label_b,
            int(previous.get("group_a_cells", -1)) == len(cells_a),
            int(previous.get("group_b_cells", -1)) == len(cells_b),
        ])
        if previous.get("status") == "complete" and identity_matches and validate_dmr(output_file):
            previous["reused"] = True
            return previous
        if previous.get("status") == "ineligible" and identity_matches:
            previous["reused"] = True
            return previous

    if len(cells_a) < args.min_cells or len(cells_b) < args.min_cells:
        record = {
            "sample_id": sample,
            "comparison": pair,
            "cell_type_a": label_a,
            "cell_type_b": label_b,
            "group_a_cells": len(cells_a),
            "group_b_cells": len(cells_b),
            "status": "ineligible",
            "dmr_rows": 0,
        }
        status_file.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return record

    with group_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows((cell, "group_A") for cell in cells_a)
        writer.writerows((cell, "group_B") for cell in cells_b)

    command = [
        str(args.methscan), "diff",
        "--threads", str(args.threads),
        "--min-cells", str(args.min_cells),
        "--bandwidth", str(args.bandwidth),
        "--stepsize", str(args.stepsize),
        "--threshold", str(args.threshold),
    ]
    if args.bridge_gaps:
        command.extend(["--bridge-gaps", str(args.bridge_gaps)])
    command.extend([str(args.data_dir), str(group_file), str(output_file)])

    with log_file.open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode != 0 or not validate_dmr(output_file):
        status = (
            "known_fdr_zero_failure"
            if completed.returncode != 0 and is_known_fdr_zero_failure(log_file)
            else "failed"
        )
        record = {
            "sample_id": sample,
            "comparison": pair,
            "cell_type_a": label_a,
            "cell_type_b": label_b,
            "group_a_cells": len(cells_a),
            "group_b_cells": len(cells_b),
            "status": status,
            "return_code": completed.returncode,
            "dmr_rows": 0,
            "group_file": str(group_file),
            "dmr_file": str(output_file),
            "log": str(log_file),
        }
        status_file.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return record

    dmr_rows = sum(1 for line in output_file.open() if line.strip())
    record = {
        "sample_id": sample,
        "comparison": pair,
        "cell_type_a": label_a,
        "cell_type_b": label_b,
        "group_a_cells": len(cells_a),
        "group_b_cells": len(cells_b),
        "status": "complete",
        "return_code": 0,
        "dmr_rows": dmr_rows,
        "group_file": str(group_file),
        "dmr_file": str(output_file),
        "log": str(log_file),
    }
    status_file.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main():
    args = parse_args()
    if not args.methscan.is_file() or not os.access(args.methscan, os.X_OK):
        raise FileNotFoundError(args.methscan)
    if not args.data_dir.is_dir():
        raise NotADirectoryError(args.data_dir)
    if args.min_cells < 1 or args.threads < 1 or args.jobs < 1:
        raise ValueError("min-cells, threads and jobs must be positive")
    allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", args.jobs * args.threads))
    if args.jobs * args.threads > allocated_cpus:
        raise ValueError(
            "jobs * threads exceeds SLURM_CPUS_PER_TASK: %d * %d > %d"
            % (args.jobs, args.threads, allocated_cpus)
        )
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.resume:
        raise FileExistsError("Meth-diff output is not empty: %s" % args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_data_dir = args.data_dir.resolve()
    args.data_dir = build_primary_data_view(source_data_dir, args.output_dir)

    metadata = read_metadata(args.cell_metadata)
    filtered_ids = read_filtered_ids(args.expected_cell_ids)
    grouped = {}
    for cell_id in filtered_ids:
        values = metadata.get(cell_id)
        if values:
            sample, cell_type = values
            if args.pool_samples_as:
                sample = args.pool_samples_as
            grouped.setdefault(sample, {}).setdefault(cell_type, []).append(cell_id)
    if not grouped or any(len(types) < 2 for types in grouped.values()):
        raise ValueError("Every sample must contain at least two eligible cell types")

    for sample, type_cells in grouped.items():
        safe_labels = [sanitize(label) for label in type_cells]
        if len(safe_labels) != len(set(safe_labels)):
            raise ValueError("Sanitized cell-type labels collide in sample %s" % sample)

    comparisons = [
        (sample, type_a, type_b, type_cells[type_a], type_cells[type_b])
        for sample, type_cells in sorted(grouped.items())
        for type_a, type_b in itertools.combinations(sorted(type_cells), 2)
    ]
    tasks = [(args, comparison) for comparison in comparisons]
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_comparison, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
            print("Completed %d/%d comparisons" % (len(results), len(tasks)), flush=True)
    results.sort(key=lambda row: (row["sample_id"], row["comparison"]))

    fields = [
        "sample_id", "comparison", "cell_type_a", "cell_type_b", "group_a_cells", "group_b_cells",
        "status", "return_code", "dmr_rows", "dmr_file", "group_file", "log",
    ]
    with (args.output_dir / "pairwise_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    summary = {
        "status": (
            "complete"
            if all(row["status"] in {"complete", "ineligible"} for row in results)
            else "complete_with_fallback_needed"
            if all(row["status"] in {"complete", "ineligible", "known_fdr_zero_failure"} for row in results)
            else "failed"
        ),
        "data_dir": str(args.data_dir.resolve()),
        "source_data_dir": str(source_data_dir),
        "chromosome_scope": "chr1-22,X,Y",
        "cell_metadata": str(args.cell_metadata.resolve()),
        "filtered_cells": len(filtered_ids),
        "eligible_cells": sum(len(cells) for types in grouped.values() for cells in types.values()),
        "samples": {
            sample: {cell_type: len(cells) for cell_type, cells in sorted(types.items())}
            for sample, types in sorted(grouped.items())
        },
        "comparisons": len(results),
        "completed_comparisons": sum(row["status"] == "complete" for row in results),
        "ineligible_comparisons": sum(row["status"] == "ineligible" for row in results),
        "failed_comparisons": sum(row["status"] == "failed" for row in results),
        "fallback_needed_comparisons": sum(
            row["status"] == "known_fdr_zero_failure" for row in results
        ),
        "min_cells": args.min_cells,
        "threads": args.threads,
        "jobs": args.jobs,
        "bandwidth": args.bandwidth,
        "stepsize": args.stepsize,
        "threshold": args.threshold,
        "sample_mode": "pooled" if args.pool_samples_as else "sample_local",
        "pooled_sample_label": args.pool_samples_as,
    }
    (args.output_dir / "pairwise_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["status"] == "failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
