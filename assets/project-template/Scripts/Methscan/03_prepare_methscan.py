#!/usr/bin/env python3
"""Run `methscan prepare --input-format bismark` from a coverage manifest."""

import argparse
import csv
import json
import subprocess
from pathlib import Path

def read_ids(path):
    with path.open() as handle:
        return [line.strip() for line in handle if line.strip()]

def read_stats_ids(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "cell_name" not in rows[0]:
        raise ValueError("MethSCAn cell_stats.csv must contain cell_name")
    return [row["cell_name"] for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methscan", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunksize", type=int, default=10000000)
    parser.add_argument("--round-sites", action="store_true")
    args = parser.parse_args()

    if not args.methscan.is_file():
        raise FileNotFoundError(args.methscan)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Prepare output is not empty: {args.output_dir}")

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "cov_path" not in rows[0]:
        raise ValueError("Coverage manifest has no input rows or cov_path column")

    paths = [row["cov_path"] for row in rows]
    missing = [path for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing coverage input, first: {missing[0]}")

    command = [
        str(args.methscan), "prepare", "--input-format", "bismark",
        "--chunksize", str(args.chunksize),
    ]
    if args.round_sites:
        command.append("--round-sites")
    command.extend(paths)
    command.append(str(args.output_dir))

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    provenance = {
        "command": command,
        "input_count": len(paths),
        "manifest": str(args.manifest.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "round_sites": args.round_sites,
        "input_format": "bismark",
    }
    with (args.output_dir.parent / "prepare_command.json").open("w") as handle:
        json.dump(provenance, handle, indent=2, sort_keys=True)
        handle.write("\n")
    subprocess.run(command, check=True)

    header_path = args.output_dir / "column_header.txt"
    stats_path = args.output_dir / "cell_stats.csv"
    if not header_path.is_file() or not stats_path.is_file():
        raise FileNotFoundError("methscan prepare did not produce column_header.txt/cell_stats.csv")
    expected_ids = sorted(row["cell_id"] for row in rows)
    header_ids = read_ids(header_path)
    stats_ids = read_stats_ids(stats_path)
    errors = []
    if header_ids != expected_ids:
        errors.append("column_header.txt != sorted coverage manifest cell_id")
    if stats_ids != expected_ids:
        errors.append("cell_stats.csv != sorted coverage manifest cell_id")
    if header_ids != stats_ids:
        errors.append("column_header.txt != cell_stats.csv")
    check = {
        "status": "pass" if not errors else "fail",
        "manifest_cells": len(expected_ids),
        "prepared_cells": len(header_ids),
        "canonical_cell_id": "<sample_id>_<17bp_barcode>",
        "errors": errors,
    }
    with (args.output_dir / "cell_id_check.json").open("w") as handle:
        json.dump(check, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if errors:
        raise ValueError("MethSCAn prepare cell-ID validation failed: %s" % "; ".join(errors))
    print(json.dumps(check, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
