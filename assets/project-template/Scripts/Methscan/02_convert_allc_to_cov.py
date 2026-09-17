#!/usr/bin/env python3
"""Convert Scanpy-selected ALLC files to CpG-only Bismark coverage (.cov.gz)."""

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "16")),
    )
    parser.add_argument("--compresslevel", type=int, default=1)
    return parser.parse_args()


def convert_one(task):
    row, cov_dir, compresslevel = task
    src = Path(row["source_path"])
    cell_id = row["cell_id"]
    out = cov_dir / f"{cell_id}.cov.gz"

    covered_cpg = 0
    methylated_reads = 0
    total_reads = 0
    methscan_n_obs = 0
    methscan_n_meth = 0

    with gzip.open(src, "rt") as inp, gzip.open(out, "wt", compresslevel=compresslevel) as dst:
        for line_number, line in enumerate(inp, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{src}:{line_number}: expected at least 6 ALLC columns")
            chrom, pos_s, _strand, context, mc_s, cov_s = fields[:6]
            # ALLCools contexts such as CGN are CpG; CH contexts are excluded.
            if not context.upper().startswith("CG"):
                continue
            try:
                pos = int(pos_s)
                mc = int(mc_s)
                cov = int(cov_s)
            except ValueError as exc:
                # Some ALLC variants contain a literal header. Skip it only when
                # it is clearly non-numeric; malformed data rows still fail below.
                if line_number == 1:
                    continue
                raise ValueError(f"{src}:{line_number}: invalid integer field") from exc
            if pos < 1 or mc < 0 or cov <= 0 or mc > cov:
                raise ValueError(f"{src}:{line_number}: invalid CpG coordinate/count")

            unmeth = cov - mc
            pct = 100.0 * mc / cov
            dst.write(f"{chrom}\t{pos}\t{pos}\t{pct:.6f}\t{mc}\t{unmeth}\n")

            covered_cpg += 1
            methylated_reads += mc
            total_reads += cov

            # MethSCAn prepare (round_sites=False) ignores mixed calls where both
            # methylated and unmethylated read counts are non-zero.
            if mc == 0 or unmeth == 0:
                methscan_n_obs += 1
                if mc > 0:
                    methscan_n_meth += 1

    if covered_cpg == 0:
        raise ValueError(f"No covered CpG records found in {src}")

    return {
        "cell_id": cell_id,
        "sample_id": row["sample_id"],
        "barcode": row["barcode"],
        "allc_path": str(src.resolve()),
        "cov_path": str(out.resolve()),
        "covered_CpG": covered_cpg,
        "methylated_CpG_reads": methylated_reads,
        "total_CpG_reads": total_reads,
        "mCG": methylated_reads / total_reads if total_reads else float("nan"),
        "methscan_n_obs_expected": methscan_n_obs,
        "methscan_n_meth_expected": methscan_n_meth,
        "methscan_global_meth_frac_expected": (
            methscan_n_meth / methscan_n_obs if methscan_n_obs else float("nan")
        ),
    }


def main():
    args = parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    if not 0 <= args.compresslevel <= 9:
        raise ValueError("compresslevel must be between 0 and 9")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Coverage output is not empty: {args.output_dir}")

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"cell_id", "sample_id", "barcode", "source_path"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Manifest must contain cell_id, sample_id, barcode and source_path")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cov_dir = args.output_dir / "cov"
    cov_dir.mkdir()
    tasks = [(row, cov_dir, args.compresslevel) for row in rows]

    results = []
    with cf.ProcessPoolExecutor(max_workers=args.workers) as executor:
        for i, result in enumerate(executor.map(convert_one, tasks), start=1):
            results.append(result)
            if i == 1 or i % 100 == 0 or i == len(tasks):
                print(f"Converted {i}/{len(tasks)} cells ({100.0 * i / len(tasks):.2f}%)", flush=True)

    manifest_fields = list(rows[0].keys()) + ["cov_path"]
    result_by_cell = {row["cell_id"]: row for row in results}
    with (args.output_dir / "input_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row, cov_path=result_by_cell[row["cell_id"]]["cov_path"]))

    qc_fields = [
        "cell_id", "sample_id", "barcode", "covered_CpG",
        "methylated_CpG_reads", "total_CpG_reads", "mCG",
        "methscan_n_obs_expected", "methscan_n_meth_expected",
        "methscan_global_meth_frac_expected", "allc_path", "cov_path",
    ]
    with (args.output_dir / "cell_cpg_qc.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    summary = {
        "status": "pass",
        "input_cells": len(rows),
        "converted_cells": len(results),
        "workers": args.workers,
        "compresslevel": args.compresslevel,
        "format": "CpG-only Bismark coverage gzip",
        "context_rule": "ALLC context startswith CG",
        "mCG_definition": "sum(methylated CpG reads) / sum(total CpG reads)",
    }
    with (args.output_dir / "conversion_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
