#!/usr/bin/env python3
"""Select all unique pooled hypo-DMRs and merge overlapping DMR intervals."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path


def read_blacklist(path: Path, chromosomes: set[str]) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            if chrom in chromosomes:
                result.setdefault(chrom, []).append((int(start), int(end)))
    for values in result.values():
        values.sort()
    return result


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
    parser.add_argument("--pairwise-summary", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--blacklist-fraction", type=float, default=0.2)
    parser.add_argument("--raw-p", type=float, default=0.01)
    parser.add_argument("--min-abs-diff", type=float, default=0.25)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sort-threads", type=int, default=8)
    args = parser.parse_args()
    if not 0 < args.raw_p <= 1 or not 0 <= args.min_abs_diff <= 1:
        raise ValueError("Invalid raw-p or methylation-difference cutoff")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    with args.chrom_sizes.open() as handle:
        chromosome_order = [line.split("\t", 1)[0].strip() for line in handle if line.strip() and not line.startswith("#")]
    if not chromosome_order or len(chromosome_order) != len(set(chromosome_order)):
        raise ValueError("chrom-sizes must contain unique chromosomes")
    ranks = {chrom: index for index, chrom in enumerate(chromosome_order)}
    blacklist = read_blacklist(args.blacklist, set(ranks))
    with args.pairwise_summary.open(newline="") as handle:
        comparisons = list(csv.DictReader(handle, delimiter="\t"))
    if not comparisons or any(row["status"] != "complete" for row in comparisons):
        raise RuntimeError("All pooled comparisons must be complete")

    unsorted = args.output_dir / "qualifying.unsorted.tsv"
    sorted_path = args.output_dir / "qualifying.sorted.tsv"
    input_rows = qualifying_rows = removed_blacklist = 0
    with unsorted.open("w") as out:
        for comparison in comparisons:
            path = Path(comparison["dmr_file"])
            if not path.is_absolute():
                path = Path.cwd() / path
            with path.open() as handle:
                for line_number, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    input_rows += 1
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 12:
                        raise ValueError(f"{path}:{line_number}: expected 12 fields")
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    if chrom not in ranks:
                        continue
                    meth_a, meth_b, raw_p = float(fields[7]), float(fields[8]), float(fields[10])
                    abs_diff = abs(meth_a - meth_b)
                    if raw_p >= args.raw_p or abs_diff < args.min_abs_diff:
                        continue
                    if blacklisted(chrom, start, end, blacklist, args.blacklist_fraction):
                        removed_blacklist += 1
                        continue
                    if fields[9] == "group_A":
                        hypo = comparison["cell_type_a"]
                    elif fields[9] == "group_B":
                        hypo = comparison["cell_type_b"]
                    else:
                        raise ValueError(f"{path}:{line_number}: invalid low group {fields[9]!r}")
                    out.write(
                        f"{ranks[chrom]}\t{chrom}\t{start}\t{end}\t{hypo}\t"
                        f"{abs_diff:.17g}\t{raw_p:.17g}\t{comparison['comparison']}\n"
                    )
                    qualifying_rows += 1

    env = os.environ.copy(); env["LC_ALL"] = "C"
    with sorted_path.open("w") as out:
        subprocess.run(
            ["sort", "--parallel", str(args.sort_threads), "-S", "50%",
             "-k1,1n", "-k3,3n", "-k4,4n", "-k5,5", str(unsorted)],
            check=True, stdout=out, env=env,
        )

    unique_gz = args.output_dir / "all_unique_hypo_DMRs.tsv.gz"
    merged_bed = args.output_dir / "all_unique_hypo_DMRs.merged.bed"
    merged_annotation = args.output_dir / "all_unique_hypo_DMRs.merged.annotation.tsv"
    unique_rows = 0
    current_key = None
    best = None
    merged = None
    merged_rows = []

    def consume(record):
        nonlocal merged
        _rank, chrom, start, end, hypo, abs_diff, raw_p, _comparison = record
        start, end, abs_diff, raw_p = int(start), int(end), float(abs_diff), float(raw_p)
        if merged is not None and chrom == merged[0] and start < merged[2]:
            merged[2] = max(merged[2], end)
            merged[3].add(hypo); merged[4] += 1
            merged[5] = max(merged[5], abs_diff); merged[6] = min(merged[6], raw_p)
        else:
            if merged is not None:
                merged_rows.append(merged)
            merged = [chrom, start, end, {hypo}, 1, abs_diff, raw_p]

    with sorted_path.open() as inp, gzip.open(unique_gz, "wt") as unique:
        unique.write("chrom\tstart\tend\thypo_cell_type\tabs_meth_diff\traw_p\tcomparison\n")
        for line in inp:
            record = line.rstrip("\n").split("\t")
            key = tuple(record[1:5])
            if key != current_key:
                if best is not None:
                    unique.write("\t".join([best[1], best[2], best[3], best[4], best[5], best[6], best[7]]) + "\n")
                    consume(best); unique_rows += 1
                current_key, best = key, record
            elif (-float(record[5]), float(record[6]), record[7]) < (-float(best[5]), float(best[6]), best[7]):
                best = record
        if best is not None:
            unique.write("\t".join([best[1], best[2], best[3], best[4], best[5], best[6], best[7]]) + "\n")
            consume(best); unique_rows += 1
    if merged is not None:
        merged_rows.append(merged)

    with merged_bed.open("w") as bed, merged_annotation.open("w", newline="") as annotation:
        writer = csv.writer(annotation, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "dmr_id", "hypo_cell_types",
                         "source_unique_dmrs", "max_abs_meth_diff", "min_raw_p"])
        for index, row in enumerate(merged_rows, 1):
            chrom, start, end, types, count, max_diff, min_p = row
            dmr_id = f"PDMR_{index:06d}"
            bed.write(f"{chrom}\t{start}\t{end}\t{max_diff:.17g}\t0\t0\t{dmr_id}\n")
            writer.writerow([chrom, start, end, dmr_id, ",".join(sorted(types)), count, max_diff, min_p])

    unsorted.unlink(); sorted_path.unlink()
    lengths = [row[2] - row[1] for row in merged_rows]
    summary = {
        "status": "complete", "pairwise_comparisons": len(comparisons),
        "input_dmr_rows": input_rows, "qualifying_rows_before_exact_dedup": qualifying_rows,
        "blacklist_removed_rows": removed_blacklist, "unique_hypo_dmrs": unique_rows,
        "merged_nonoverlapping_dmrs": len(merged_rows),
        "merged_total_bp": sum(lengths), "merged_max_length": max(lengths, default=0),
        "raw_p_strictly_less_than": args.raw_p, "min_abs_meth_diff": args.min_abs_diff,
        "top_n_truncation": None, "merged_bed": str(merged_bed.resolve()),
        "annotation": str(merged_annotation.resolve()),
    }
    (args.output_dir / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
