#!/usr/bin/env python3
"""Prepare a canonical, blacklist-filtered VMR BED and a validated ALLC table."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path


COV_SUFFIX = "_allc.gz.cov"


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def read_chrom_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    order, sizes = [], {}
    with path.open() as handle:
        for line in handle:
            chrom, size_s = line.rstrip().split("\t")[:2]
            if chrom in sizes:
                raise ValueError(f"Duplicate chromosome in {path}: {chrom}")
            order.append(chrom)
            sizes[chrom] = int(size_s)
    return order, sizes


def read_blacklist(path: Path, chromosomes: set[str]) -> dict[str, list[tuple[int, int]]]:
    regions = {chrom: [] for chrom in chromosomes}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom in regions:
                regions[chrom].append((start, end))
    for values in regions.values():
        values.sort()
    return regions


def blacklisted(chrom: str, start: int, end: int, regions, fraction: float) -> bool:
    length = end - start
    for left, right in regions.get(chrom, ()):
        if left >= end:
            break
        if right <= start:
            continue
        if (min(end, right) - max(start, left)) / length >= fraction:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-bed", type=Path, default=Path(os.environ["VMR_SOURCE_BED"]))
    parser.add_argument("--output-bed", type=Path, default=Path(os.environ["VMR_FILTERED_BED"]))
    parser.add_argument("--chrom-sizes", type=Path, default=Path(os.environ["VMR_CHROM_SIZES"]))
    parser.add_argument("--blacklist", type=Path, default=Path(os.environ["VMR_BLACKLIST"]))
    parser.add_argument("--blacklist-md5", default=os.environ["VMR_BLACKLIST_MD5"])
    parser.add_argument("--blacklist-fraction", type=float, default=float(os.environ["VMR_BLACKLIST_FRACTION"]))
    parser.add_argument("--cov-dir", type=Path, default=Path(os.environ["VMR_COV_DIR"]))
    parser.add_argument("--allc-dir", type=Path, default=Path(os.environ["VMR_EXISTING_ALLC_DIR"]))
    parser.add_argument(
        "--input-manifest", type=Path, default=Path(os.environ.get("VMR_INPUT_MANIFEST", "")),
        help="MethSCAn selected-ALLC manifest; takes precedence over the legacy coverage fallback",
    )
    parser.add_argument(
        "--filtered-cell-ids", type=Path,
        default=(Path(os.environ["VMR_FILTERED_CELL_IDS"])
                 if os.environ.get("VMR_FILTERED_CELL_IDS") else None),
        help="MethSCAn filter column_header.txt; restricts the manifest to retained cells",
    )
    parser.add_argument("--allc-table", type=Path, default=Path(os.environ["VMR_ALLC_TABLE"]))
    parser.add_argument("--expected-cells", type=int, default=int(os.environ["VMR_EXPECTED_CELLS"]))
    args = parser.parse_args()
    if not 0 < args.blacklist_fraction <= 1:
        raise ValueError("blacklist fraction must be in (0, 1]")
    for path in (args.source_bed, args.chrom_sizes, args.blacklist):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    observed_md5 = digest(args.blacklist, "md5")
    if args.blacklist_md5 and observed_md5.lower() != args.blacklist_md5.lower():
        raise ValueError(f"Blacklist MD5 mismatch: {observed_md5}")

    chrom_order, chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    order_index = {chrom: index for index, chrom in enumerate(chrom_order)}
    blacklist = read_blacklist(args.blacklist, set(chrom_sizes))
    source_count = noncanonical = removed_blacklist = 0
    retained: list[tuple[str, int, int, float, int, int]] = []
    with open_text(args.source_bed) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 6:
                raise ValueError(f"{args.source_bed}:{line_number}: expected MethSCAn 6-column BED")
            source_count += 1
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in chrom_sizes:
                noncanonical += 1
                continue
            if start < 0 or end <= start or end > chrom_sizes[chrom]:
                raise ValueError(f"{args.source_bed}:{line_number}: invalid interval")
            if blacklisted(chrom, start, end, blacklist, args.blacklist_fraction):
                removed_blacklist += 1
                continue
            retained.append((chrom, start, end, float(fields[3]), int(fields[4]), int(fields[5])))
    retained.sort(key=lambda row: (order_index[row[0]], row[1], row[2]))
    for previous, current in zip(retained, retained[1:]):
        if previous[0] == current[0] and current[1] < previous[2]:
            raise ValueError(f"Overlapping VMRs are not supported: {previous} and {current}")

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    bed_text = "".join(
        f"{chrom}\t{start}\t{end}\t{peak_var:.17g}\t{n_cpg}\t{n_obs_cells}\tVMR_{index:06d}\n"
        for index, (chrom, start, end, peak_var, n_cpg, n_obs_cells) in enumerate(retained, start=1)
    )
    temporary_bed = args.output_bed.with_suffix(".tmp.bed")
    temporary_bed.write_text(bed_text)
    temporary_bed.replace(args.output_bed)

    selected: list[tuple[str, Path, Path]] = []
    if args.input_manifest.is_file():
        filtered_ids = None
        if args.filtered_cell_ids is not None:
            if not args.filtered_cell_ids.is_file():
                raise FileNotFoundError(args.filtered_cell_ids)
            filtered_ids = {line.strip() for line in args.filtered_cell_ids.open() if line.strip()}
            if not filtered_ids:
                raise ValueError(f"Filtered cell list is empty: {args.filtered_cell_ids}")
        with args.input_manifest.open(newline="") as handle:
            manifest = csv.DictReader(handle, delimiter="\t")
            required = {"cell_id", "source_path", "source_index"}
            if not manifest.fieldnames or not required.issubset(manifest.fieldnames):
                raise ValueError(f"{args.input_manifest} requires columns {sorted(required)}")
            for row in manifest:
                cell, allc, index = row["cell_id"], Path(row["source_path"]), Path(row["source_index"])
                if filtered_ids is not None and cell not in filtered_ids:
                    continue
                if not cell or not allc.is_file() or allc.stat().st_size == 0 or not index.is_file() or index.stat().st_size == 0:
                    raise FileNotFoundError(f"Invalid selected ALLC for {cell}: {allc}")
                selected.append((cell, allc.resolve(), index.resolve()))
        if not selected or len({cell for cell, _allc, _index in selected}) != len(selected):
            raise ValueError(f"{args.input_manifest} has no cells or duplicate cell IDs")
        if filtered_ids is not None and {cell for cell, _allc, _index in selected} != filtered_ids:
            missing = sorted(filtered_ids - {cell for cell, _allc, _index in selected})
            raise ValueError(f"Filtered cells missing from MethSCAn manifest (first: {missing[0] if missing else 'unknown'})")
        source_kind = f"MethSCAn manifest: {args.input_manifest.resolve()}"
    else:
        cov_files = sorted(args.cov_dir.glob(f"*{COV_SUFFIX}"))
        if args.expected_cells and len(cov_files) != args.expected_cells:
            raise RuntimeError(f"Expected {args.expected_cells} coverage files, found {len(cov_files)}")
        for cov in cov_files:
            cell = cov.name[: -len(COV_SUFFIX)]
            allc = args.allc_dir / f"{cell}.allc.tsv.gz"
            index = Path(f"{allc}.tbi")
            if not allc.is_file() or allc.stat().st_size == 0 or not index.is_file() or index.stat().st_size == 0:
                raise FileNotFoundError(f"Missing ALLC or tabix index for {cell}: {allc}")
            selected.append((cell, allc.resolve(), index.resolve()))
        if not selected:
            raise RuntimeError(f"No coverage files in {args.cov_dir}")
        source_kind = f"legacy coverage directory: {args.cov_dir.resolve()}"
    rows, manifest_rows = [], []
    for cell, allc, index in selected:
        rows.append(f"{cell}\t{allc}\n")
        manifest_rows.append(
            f"{cell}\t{allc}\t{allc.stat().st_size}\t{allc.stat().st_mtime_ns}\t"
            f"{index.stat().st_size}\t{index.stat().st_mtime_ns}\n"
        )
    args.allc_table.parent.mkdir(parents=True, exist_ok=True)
    temporary_table = args.allc_table.with_suffix(".tmp.tsv")
    temporary_table.write_text("".join(rows))
    temporary_table.replace(args.allc_table)
    manifest = args.allc_table.parent / "input_allc_manifest.tsv"
    manifest.write_text(
        "cell_id\tpath\tsize\tmtime_ns\tindex_size\tindex_mtime_ns\n" + "".join(manifest_rows)
    )
    summary = {
        "source_bed": str(args.source_bed.resolve()), "source_vmr_count": source_count,
        "noncanonical_removed": noncanonical, "blacklist_removed": removed_blacklist,
        "retained_vmr_count": len(retained), "blacklist_fraction": args.blacklist_fraction,
        "blacklist_md5": observed_md5, "cells": len(rows), "allc_source": source_kind,
        "filtered_bed": str(args.output_bed.resolve()), "allc_table": str(args.allc_table.resolve()),
    }
    (args.output_bed.parent / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
