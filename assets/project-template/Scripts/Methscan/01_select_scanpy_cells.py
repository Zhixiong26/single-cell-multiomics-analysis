#!/usr/bin/env python3
"""Discover ALLCs from explicit sample rows and retain approved annotated cells."""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-tsv", type=Path)
    parser.add_argument("--allc-source", type=Path)
    parser.add_argument("--scanpy-annotation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample", action="append", dest="samples")
    parser.add_argument("--exclude-cell-type", default="NA")
    parser.add_argument("--annotation-approved", action="store_true")
    return parser.parse_args()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def derive_cell_name(path, row):
    regex = (row.get("allc_cell_id_regex") or "").strip()
    replacement = row.get("allc_cell_id_replacement") or ""
    if regex:
        match = re.search(regex, path.name)
        if not match:
            raise ValueError("ALLC filename does not match configured cell-ID regex: %s" % path)
        value = match.expand(replacement) if replacement else (match.groupdict().get("cell_id") or match.group(0))
    else:
        value = path.name
        for suffix in (".allc.tsv.gz", "_allc.gz", ".allc.gz"):
            if value.endswith(suffix):
                value = value[:-len(suffix)]
                break
    prefix = (row.get("cell_id_prefix") or row["sample_id"]).strip()
    cell_id = value if not prefix or value.startswith(prefix + "_") else prefix + "_" + value
    return value, cell_id


def discover_allcs(sample_rows):
    rows = []
    seen = set()
    for sample in sample_rows:
        source = Path(sample["allc_root"]).resolve()
        pattern = (sample.get("allc_glob") or "**/*.allc.tsv.gz").strip()
        paths = sorted(path for path in source.glob(pattern) if path.is_file())
        if not paths:
            raise ValueError("No ALLC files found for sample %s with %s" % (sample["sample_id"], pattern))
        for path in paths:
            barcode, cell_id = derive_cell_name(path, sample)
            if cell_id in seen:
                raise ValueError("Duplicate ALLC cell_id discovered: %s" % cell_id)
            seen.add(cell_id)
            index = Path(str(path) + ".tbi")
            if not index.is_file():
                raise FileNotFoundError("Missing ALLC index: %s" % index)
            rows.append({"sample_id": sample["sample_id"], "barcode": barcode,
                         "cell_id": cell_id, "source_path": str(path.resolve()),
                         "source_index": str(index.resolve())})
    return rows


def main():
    args = parse_args()
    if not args.annotation_approved:
        raise ValueError("MethSCAn cell-type routes require approved annotation")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Scanpy selection output is not empty: %s" % args.output_dir)

    if args.samples_tsv:
        sample_rows = [row for row in read_tsv(args.samples_tsv)
                       if row.get("include", "").lower() in {"1", "true", "yes"} and row.get("allc_root")]
    elif args.allc_source and args.samples:
        sample_rows = [{"sample_id": sample, "allc_root": str(args.allc_source),
                        "allc_glob": "**/*.gz", "cell_id_prefix": sample}
                       for sample in args.samples]
    else:
        raise ValueError("provide --samples-tsv or legacy --allc-source/--sample arguments")
    requested_samples = {row["sample_id"] for row in sample_rows}
    allc_rows = discover_allcs(sample_rows)
    if not allc_rows:
        raise ValueError("No matching per-cell ALLC files found")

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
        elif cell_type in {args.exclude_cell_type, "NA", "Unassigned", "requires_review"}:
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
