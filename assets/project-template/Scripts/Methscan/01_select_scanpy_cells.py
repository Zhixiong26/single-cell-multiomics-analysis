#!/usr/bin/env python3
"""Discover ALLC files and keep cells present in the non-NA Scanpy whitelist."""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


ALLC_NAME_RE = re.compile(
    r"^(?:(?P<sample>[A-Za-z0-9]+)_)?"
    r"(?P<barcode>[ACGT]{17})"
    r"(?:_allc|\.allc(?:\.tsv)?)\.gz$"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allc-source", type=Path, required=True)
    parser.add_argument("--scanpy-annotation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample", action="append", dest="samples", required=True)
    parser.add_argument("--exclude-cell-type", default="NA")
    return parser.parse_args()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def infer_sample(path, source, requested_samples, filename_sample=None):
    if filename_sample in requested_samples:
        return filename_sample

    relative = path.relative_to(source)
    hits = set()
    for part in relative.parts[:-1]:
        tokens = re.split(r"[^A-Za-z0-9]+", part)
        hits.update(sample for sample in requested_samples if sample in tokens)

    if len(hits) == 1:
        return next(iter(hits))
    if len(hits) > 1:
        raise ValueError("Ambiguous sample for ALLC path: %s" % path)
    return None


def discover_allcs(source, requested_samples):
    rows = []
    seen = set()

    for path in sorted(source.rglob("*.gz")):
        match = ALLC_NAME_RE.match(path.name)
        if not match:
            continue

        sample = infer_sample(
            path,
            source,
            requested_samples,
            filename_sample=match.group("sample"),
        )
        if sample is None:
            continue

        barcode = match.group("barcode")
        cell_id = "%s_%s" % (sample, barcode)
        if cell_id in seen:
            raise ValueError("Duplicate ALLC cell_id discovered: %s" % cell_id)
        seen.add(cell_id)

        index = Path(str(path) + ".tbi")
        rows.append({
            "sample_id": sample,
            "barcode": barcode,
            "cell_id": cell_id,
            "source_path": str(path.resolve()),
            "source_index": str(index.resolve()) if index.is_file() else "NA",
        })

    return rows


def main():
    args = parse_args()
    source = args.allc_source.resolve()
    requested_samples = set(args.samples)

    if not source.is_dir():
        raise NotADirectoryError(source)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Scanpy selection output is not empty: %s" % args.output_dir)

    allc_rows = discover_allcs(source, requested_samples)
    if not allc_rows:
        raise ValueError("No matching per-cell ALLC files found below %s" % source)

    discovered_by_sample = Counter(row["sample_id"] for row in allc_rows)
    missing_samples = sorted(requested_samples.difference(discovered_by_sample))
    if missing_samples:
        raise ValueError(
            "No ALLC files found for required sample(s): %s" % ", ".join(missing_samples)
        )

    annotation_rows = read_tsv(args.scanpy_annotation)
    required = {"cell_id", "cell_type"}
    if not annotation_rows or not required.issubset(annotation_rows[0]):
        raise ValueError("Scanpy annotation must contain cell_id and cell_type")
    if len({row["cell_id"] for row in annotation_rows}) != len(annotation_rows):
        raise ValueError("Scanpy annotation cell_id values must be unique")
    annotation = {row["cell_id"]: row for row in annotation_rows}

    selected, excluded = [], []
    for row in allc_rows:
        annotation_row = annotation.get(row["cell_id"])
        if annotation_row is None:
            excluded.append({
                "cell_id": row["cell_id"],
                "sample_id": row["sample_id"],
                "source_path": row["source_path"],
                "reason": "not_in_scanpy",
            })
            continue

        cell_type = (annotation_row.get("cell_type") or "").strip()
        if not cell_type:
            reason = "missing_scanpy_cell_type"
        elif cell_type == args.exclude_cell_type:
            reason = "excluded_scanpy_cell_type_%s" % args.exclude_cell_type
        else:
            reason = None

        if reason:
            excluded.append({
                "cell_id": row["cell_id"],
                "sample_id": row["sample_id"],
                "source_path": row["source_path"],
                "reason": reason,
            })
            continue

        selected.append(dict(row, **{
            "rna_cohort": annotation_row.get("cohort", ""),
            "rna_sample": annotation_row.get("sample", ""),
            "rna_leiden": annotation_row.get("leiden", ""),
            "rna_cell_type": cell_type,
            "rna_annotation_available": "True",
        }))

    if not selected:
        raise ValueError("No ALLC cells remain after Scanpy non-NA selection")

    selected_by_sample = Counter(row["sample_id"] for row in selected)
    missing_selected_samples = sorted(requested_samples.difference(selected_by_sample))
    if missing_selected_samples:
        raise ValueError(
            "Scanpy selection retained no cells for required sample(s): %s"
            % ", ".join(missing_selected_samples)
        )

    args.output_dir.mkdir(parents=True)

    fields = [
        "sample_id", "barcode", "cell_id", "source_path", "source_index",
        "rna_cohort", "rna_sample", "rna_leiden", "rna_cell_type",
        "rna_annotation_available",
    ]
    with (args.output_dir / "input_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    with (args.output_dir / "allc_excluded_by_scanpy.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["cell_id", "sample_id", "source_path", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(excluded)

    summary = {
        "allc_source": str(source),
        "discovered_allc_cells": len(allc_rows),
        "discovered_by_sample": dict(sorted(discovered_by_sample.items())),
        "scanpy_annotation_cells": len(annotation_rows),
        "scanpy_selected_cells": len(selected),
        "selected_by_sample": dict(sorted(selected_by_sample.items())),
        "excluded_cells": len(excluded),
        "excluded_by_reason": dict(
            sorted(Counter(row["reason"] for row in excluded).items())
        ),
        "excluded_scanpy_cell_type": args.exclude_cell_type,
    }
    with (args.output_dir / "scanpy_selection_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
