#!/usr/bin/env python3
"""Resolve known MethSCAn FDR failures, select Top200 hypo-DMRs, and plot them.

The implementation follows the referenced workflow: comparisons are sample-local;
raw p < 0.01 and absolute mean methylation difference >= 0.25 are required; at
most 200 unique intervals are retained per sample/cell type.  Single-cell values
are the arithmetic mean of unique-CpG methylation ratios within each merged DMR.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PRIMARY = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")


@dataclass(frozen=True)
class Record:
    chrom: str
    start: int
    end: int
    sample: str
    hypo: str
    other: str
    abs_diff: float
    raw_p: float
    columns: tuple[str, ...]
    source: str


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    dmr_id: str
    hypo_types: tuple[str, ...]
    source_count: int


def args_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--methscan", type=Path, required=True)
    p.add_argument("--methscan-python", type=Path, required=True)
    p.add_argument("--raw-p", type=float, default=0.01)
    p.add_argument("--min-abs-diff", type=float, default=0.25)
    p.add_argument("--top-per-cell-type", type=int, default=200)
    p.add_argument("--min-cells", type=int, default=6)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--bandwidth", type=int, default=2000)
    p.add_argument("--stepsize", type=int, default=1000)
    p.add_argument("--threshold", type=float, default=0.02)
    p.add_argument("--matrix-workers", type=int, default=32)
    p.add_argument("--zscore-min-observed-cells", type=int, default=30)
    p.add_argument("--zscore-clip", type=float, default=3.0)
    p.add_argument("--compressed-colorbar-limit", type=float, default=1.0)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def safe(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chrom_key(chrom):
    suffix = chrom[3:] if chrom.startswith("chr") else chrom
    if suffix.isdigit():
        return (0, int(suffix))
    if suffix == "X":
        return (0, 23)
    if suffix == "Y":
        return (0, 24)
    return (1, chrom)


def stage_fallback(args, comparisons):
    root = args.run_dir / "08_rawp_fallback"
    root.mkdir()
    pending = [r for r in comparisons if r["status"] == "known_fdr_zero_failure"]
    patched_root = root / "patched_python_package"
    if pending:
        result = subprocess.run(
            [str(args.methscan_python), "-c", "import os,methscan;print(os.path.dirname(methscan.__file__))"],
            check=True, capture_output=True, text=True,
        )
        source_pkg = Path(result.stdout.strip())
        shutil.copytree(source_pkg, patched_root / "methscan")
        diff_py = patched_root / "methscan" / "diff.py"
        text = diff_py.read_text()
        needle = '    adj_p_val = calc_fdr(output_final[11] == "real")\n'
        replacement = '''    # RAW_P_FALLBACK_NO_NULL_DMRS: retain raw p when null DMRs are absent.\n    is_real = output_final[11] == "real"\n    if not np.any(is_real):\n        raise RuntimeError("raw-p fallback unavailable: no real DMRs")\n    if not np.any(~is_real):\n        adj_p_val = np.full(is_real.shape, np.nan, dtype=np.float64)\n    else:\n        adj_p_val = calc_fdr(is_real)\n'''
        if text.count(needle) != 1:
            raise RuntimeError("Cannot uniquely patch MethSCAn calc_fdr call")
        diff_py.write_text(text.replace(needle, replacement))

    pairwise_summary = json.loads((args.run_dir / "07_methdiff" / "pairwise_summary.json").read_text())
    diff_data_dir = Path(pairwise_summary["data_dir"])
    results = []
    for row in comparisons:
        if row["status"] != "known_fdr_zero_failure":
            continue
        pair_dir = root / "samples" / safe(row["sample_id"]) / "comparisons" / row["comparison"]
        pair_dir.mkdir(parents=True)
        raw_output = pair_dir / "DMRs.raw.bed"
        output = pair_dir / "DMRs.bed"
        log = pair_dir / "methscan_diff.log"
        command = [
            str(args.methscan), "diff", "--threads", str(args.threads),
            "--min-cells", str(args.min_cells), "--bandwidth", str(args.bandwidth),
            "--stepsize", str(args.stepsize), "--threshold", str(args.threshold),
            str(diff_data_dir), row["group_file"], str(raw_output),
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = str(patched_root) + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        with log.open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env)
        if completed.returncode:
            raise RuntimeError(f"Raw-p fallback failed for {row['sample_id']} {row['comparison']}; see {log}")
        count = 0
        with raw_output.open() as inp, output.open("w") as out:
            for number, line in enumerate(inp, 1):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 12 or fields[11].lower() != "nan":
                    raise ValueError(f"{raw_output}:{number}: invalid fallback row")
                fields[11] = "NA"
                out.write("\t".join(fields) + "\n")
                count += 1
        raw_output.unlink()
        results.append({"sample_id": row["sample_id"], "comparison": row["comparison"],
                        "status": "complete", "dmr_rows": count, "dmr_file": str(output)})

    summary = {"status": "complete", "fallback_needed": len(pending),
               "fallback_completed": len(results), "comparisons": results}
    (root / "fallback_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    with (root / "fallback_summary.tsv").open("w", newline="") as handle:
        fields = ["sample_id", "comparison", "status", "dmr_rows", "dmr_file"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(results)
    return {(r["sample_id"], r["comparison"]): r["dmr_file"] for r in results}, summary


def parse_dmr_file(path, row, args):
    records = []
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 12:
                raise ValueError(f"{path}:{number}: expected 12 columns")
            if not PRIMARY.fullmatch(f[0]):
                continue
            start, end = int(f[1]), int(f[2])
            meth_a, meth_b, raw_p = float(f[7]), float(f[8]), float(f[10])
            if raw_p >= args.raw_p or abs(meth_a - meth_b) < args.min_abs_diff:
                continue
            if f[9] == "group_A":
                hypo, other = row["cell_type_a"], row["cell_type_b"]
            elif f[9] == "group_B":
                hypo, other = row["cell_type_b"], row["cell_type_a"]
            else:
                raise ValueError(f"{path}:{number}: invalid low_group_label {f[9]!r}")
            records.append(Record(f[0], start, end, row["sample_id"], hypo, other,
                                  abs(meth_a - meth_b), raw_p, tuple(f), str(path)))
    return records


def merge_regions(sample, selected):
    flat = [(r.chrom, r.start, r.end, r.hypo) for values in selected.values() for r in values]
    flat.sort(key=lambda x: (chrom_key(x[0]), x[1], x[2], x[3]))
    merged = []
    if not flat:
        return merged
    chrom, start, end, hypo = flat[0]
    types, count = {hypo}, 1
    for next_chrom, next_start, next_end, next_hypo in flat[1:]:
        if next_chrom == chrom and next_start < end:
            end = max(end, next_end); types.add(next_hypo); count += 1
        else:
            merged.append(Region(chrom, start, end, "", tuple(sorted(types)), count))
            chrom, start, end, types, count = next_chrom, next_start, next_end, {next_hypo}, 1
    merged.append(Region(chrom, start, end, "", tuple(sorted(types)), count))
    return [Region(r.chrom, r.start, r.end, f"{sample}__merged_hypo_DMR_{i:06d}",
                   r.hypo_types, r.source_count) for i, r in enumerate(merged, 1)]


def stage_select(args, comparisons, fallback_paths):
    root = args.run_dir / "09_top200_hypo_DMRs"
    root.mkdir()
    grouped = defaultdict(list)
    input_rows = 0
    for row in comparisons:
        if row["status"] == "ineligible":
            continue
        if row["status"] == "known_fdr_zero_failure":
            path = Path(fallback_paths[(row["sample_id"], row["comparison"])])
        elif row["status"] == "complete":
            path = Path(row["dmr_file"])
        else:
            raise RuntimeError(f"Unresolved pairwise comparison: {row}")
        records = parse_dmr_file(path, row, args)
        input_rows += int(row.get("dmr_rows") or 0)
        for record in records:
            grouped[(record.sample, record.hypo)].append(record)

    samples = sorted({r["sample_id"] for r in comparisons})
    sample_summaries = []
    for sample in samples:
        sample_dir = root / "samples" / safe(sample)
        by_type = sample_dir / "by_cell_type"
        by_type.mkdir(parents=True)
        selected = {}
        for (record_sample, cell_type), records in sorted(grouped.items()):
            if record_sample != sample:
                continue
            best = {}
            for record in records:
                key = (record.chrom, record.start, record.end)
                current = best.get(key)
                if current is None or (-record.abs_diff, record.raw_p, record.other) < (-current.abs_diff, current.raw_p, current.other):
                    best[key] = record
            ranked = sorted(best.values(), key=lambda r: (-r.abs_diff, r.raw_p, chrom_key(r.chrom), r.start, r.end))
            selected[cell_type] = ranked[:args.top_per_cell_type]
            out = by_type / f"{safe(cell_type)}__top{args.top_per_cell_type}_hypo_DMRs.tsv"
            with out.open("w", newline="") as handle:
                fields = ["chrom", "start", "end", "hypo_cell_type", "other_cell_type", "abs_meth_diff", "raw_p", "source_file"]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
                for r in selected[cell_type]:
                    writer.writerow({"chrom": r.chrom, "start": r.start, "end": r.end,
                                     "hypo_cell_type": r.hypo, "other_cell_type": r.other,
                                     "abs_meth_diff": r.abs_diff, "raw_p": r.raw_p, "source_file": r.source})
        merged = merge_regions(sample, selected)
        bed = sample_dir / "merged_top200_hypo_DMRs.bed"
        annotation = sample_dir / "merged_top200_hypo_DMRs_annotation.tsv"
        with bed.open("w") as handle:
            for r in merged: handle.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.dmr_id}\n")
        with annotation.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "start", "end", "dmr_id", "hypo_cell_types", "source_region_count"])
            for r in merged: writer.writerow([r.chrom, r.start, r.end, r.dmr_id, ",".join(r.hypo_types), r.source_count])
        sample_summaries.append({"sample_id": sample, "cell_types_with_hypo_dmrs": len(selected),
                                 "selected_unique_intervals": sum(map(len, selected.values())),
                                 "merged_dmrs": len(merged), "annotation": str(annotation)})
    summary = {"status": "complete", "raw_p_strictly_less_than": args.raw_p,
               "min_abs_methylation_difference": args.min_abs_diff,
               "top_unique_dmrs_per_cell_type": args.top_per_cell_type,
               "source_dmr_rows": input_rows, "samples": sample_summaries}
    (root / "selection_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


_INTERVALS = None
_N_REGIONS = 0


def init_worker(index, n):
    global _INTERVALS, _N_REGIONS
    _INTERVALS, _N_REGIONS = index, n


def cell_means(payload):
    cell_id, cov_path = payload
    sums = np.zeros(_N_REGIONS, dtype=np.float64)
    counts = np.zeros(_N_REGIONS, dtype=np.uint32)
    current = None
    mc = uc = 0.0
    rows = 0
    completed_chroms = set()

    def finish():
        if current is None:
            return
        chrom, pos = current
        data = _INTERVALS.get(chrom)
        if data is None:
            return
        starts, ends, indices = data
        j = bisect_right(starts, pos) - 1
        if j >= 0 and pos < ends[j]:
            total = mc + uc
            if total <= 0:
                raise ValueError(f"{cov_path}: zero coverage at {chrom}:{pos}")
            sums[indices[j]] += mc / total; counts[indices[j]] += 1

    with gzip.open(cov_path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, methylated, unmethylated = f[0], int(f[1]), float(f[4]), float(f[5])
            key = (chrom, pos)
            if current is not None and chrom != current[0]:
                completed_chroms.add(current[0])
                if chrom in completed_chroms:
                    raise ValueError(f"{cov_path}:{line_number}: chromosome {chrom} reappears")
            elif current is not None and chrom == current[0] and pos < current[1]:
                raise ValueError(f"{cov_path}:{line_number}: unsorted coordinate {chrom}:{pos}")
            if current != key:
                finish(); current, mc, uc, rows = key, methylated, unmethylated, 1
            else:
                mc += methylated; uc += unmethylated; rows += 1
    finish()
    values = np.full(_N_REGIONS, np.nan, dtype=np.float32)
    observed = counts > 0
    values[observed] = sums[observed] / counts[observed]
    return cell_id, values


def read_regions(path):
    regions = []
    for row in read_tsv(path):
        regions.append(Region(row["chrom"], int(row["start"]), int(row["end"]), row["dmr_id"],
                              tuple(filter(None, row["hypo_cell_types"].split(","))), int(row["source_region_count"])))
    return regions


def stage_matrix(args, selection):
    root = args.run_dir / "10_single_cell_DMR_matrix"
    root.mkdir()
    metadata = {r["cell_id"]: r for r in read_tsv(args.run_dir / "01_cov" / "input_manifest.tsv")}
    filtered = {line.strip() for line in (args.run_dir / "03_filtered" / "column_header.txt").open() if line.strip()}
    summaries = []
    for sample_row in selection["samples"]:
        sample = sample_row["sample_id"]
        regions = read_regions(Path(sample_row["annotation"]))
        cells = sorted((cid, r) for cid, r in metadata.items() if cid in filtered and r["sample_id"] == sample)
        if not cells:
            raise ValueError(f"{sample}: no filtered cells with COV metadata")
        sample_dir = root / "samples" / safe(sample); sample_dir.mkdir(parents=True)
        output = sample_dir / "single_cell_DMR_mean_unique_CpG_ratio.tsv.gz"
        if not regions:
            with gzip.open(output, "wt") as handle: handle.write("chrom\tstart\tend\tdmr_id\n")
            summaries.append({"sample_id": sample, "dmrs": 0, "single_cells": len(cells), "covered_values": 0, "matrix_file": str(output)})
            continue
        index = {}
        for i, r in enumerate(regions):
            starts, ends, indices = index.setdefault(r.chrom, ([], [], []))
            starts.append(r.start); ends.append(r.end); indices.append(i)
        matrix = np.full((len(regions), len(cells)), np.nan, dtype=np.float32)
        tasks = [(cid, row["cov_path"]) for cid, row in cells]
        with ProcessPoolExecutor(max_workers=min(args.matrix_workers, len(tasks)), initializer=init_worker,
                                 initargs=(index, len(regions))) as executor:
            for done, (cell_id, values) in enumerate(executor.map(cell_means, tasks), 1):
                matrix[:, done - 1] = values
                if done == 1 or done % 100 == 0 or done == len(tasks):
                    print(f"[{sample}] DMR matrix cells {done}/{len(tasks)}", flush=True)
        with gzip.open(output, "wt", newline="", compresslevel=3) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "start", "end", "dmr_id", *[c[0] for c in cells]])
            for i, r in enumerate(regions):
                writer.writerow([r.chrom, r.start, r.end, r.dmr_id, *["NA" if not np.isfinite(x) else f"{x:.8g}" for x in matrix[i]]])
        annotations = sample_dir / "cell_annotations.tsv"
        with annotations.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["cell_id", "sample_id", "cell_type", "cov_path"])
            for cid, row in cells: writer.writerow([cid, sample, row["rna_cell_type"], row["cov_path"]])
        summaries.append({"sample_id": sample, "dmrs": len(regions), "single_cells": len(cells),
                          "covered_values": int(np.isfinite(matrix).sum()), "matrix_file": str(output),
                          "cell_annotations": str(annotations), "dmr_annotations": sample_row["annotation"]})
    summary = {"status": "complete", "value_definition": "unweighted arithmetic mean of unique-CpG ratios",
               "interval_rule": "BED half-open: start <= CpG position < end", "samples": summaries}
    (root / "matrix_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def load_matrix(path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader); rows = list(reader)
    values = np.asarray([[np.nan if x == "NA" else float(x) for x in row[4:]] for row in rows], dtype=np.float32).T
    return header[4:], rows, values


def zscore(matrix, min_observed):
    observed = np.isfinite(matrix); counts = observed.sum(axis=0)
    means = np.nanmean(matrix, axis=0)
    centered = matrix - means
    sd = np.sqrt(np.nanmean(centered * centered, axis=0))
    eligible = (counts >= min_observed) & (sd > 0)
    out = np.full(matrix.shape, np.nan, dtype=np.float32)
    np.divide(centered, sd, out=out, where=observed & eligible[np.newaxis, :])
    return out, counts, eligible


def plot_one(matrix, row_types, col_types, output, title, vmin, vmax, raw, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    labels = sorted(set(row_types) | set(col_types))
    palette = plt.get_cmap("tab20", max(1, len(labels)))
    colors = {label: palette(i) for i, label in enumerate(labels)}

    def group_layout(values):
        """Return the center, label, and end boundary of contiguous groups."""
        groups = []
        start = 0
        for index in range(1, len(values) + 1):
            if index == len(values) or values[index] != values[start]:
                groups.append(((start + index - 1) / 2, values[start], index - 0.5))
                start = index
        return groups

    row_groups = group_layout(row_types)
    col_groups = group_layout(col_types)
    width = max(12, min(42, matrix.shape[1] / 70))
    height = max(10, min(42, matrix.shape[0] / 160))
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=(0.025, 1), width_ratios=(0.018, 1, 0.025), hspace=.01, wspace=.02)
    top, left, ax, cax = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])
    type_cmap = matplotlib.colors.ListedColormap([colors[x] for x in labels])
    codes = {x: i for i, x in enumerate(labels)}
    top.imshow(np.asarray([codes[x] for x in col_types])[None, :], aspect="auto", cmap=type_cmap, vmin=-.5, vmax=len(labels)-.5)
    left.imshow(np.asarray([codes[x] for x in row_types])[:, None], aspect="auto", cmap=type_cmap, vmin=-.5, vmax=len(labels)-.5)
    cmap = LinearSegmentedColormap.from_list("meth", ["#0000ff", "#ff0000"] if raw else ["#2166ac", "#f7f7f7", "#b2182b"])
    cmap.set_bad("#d9d9d9")
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    top.set_xticks([]); top.set_yticks([])
    left.set_xticks([])
    left.set_yticks([center for center, _, _ in row_groups])
    left.set_yticklabels([label for _, label, _ in row_groups], fontsize=8)
    left.tick_params(axis="y", which="both", length=0, pad=4)
    ax.set_yticks([])
    ax.set_xticks([center for center, _, _ in col_groups])
    ax.set_xticklabels([label for _, label, _ in col_groups], rotation=90,
                       ha="center", va="top", fontsize=8)
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    for _, _, boundary in row_groups[:-1]:
        ax.axhline(boundary, color="black", linewidth=0.45, alpha=0.65)
        left.axhline(boundary, color="black", linewidth=0.45, alpha=0.65)
    for _, _, boundary in col_groups[:-1]:
        ax.axvline(boundary, color="black", linewidth=0.45, alpha=0.65)
        top.axvline(boundary, color="black", linewidth=0.45, alpha=0.65)
    ax.set_title(title)
    ax.set_xlabel(
        f"All DMRs grouped by supporting hypo cell type "
        f"({matrix.shape[1]:,} exact DMR columns)", labelpad=8
    )
    ax.set_ylabel(
        f"Single cells grouped by annotated cell type ({matrix.shape[0]:,} cells)",
        labelpad=8,
    )
    ax.text(1.0, 1.005, "Gray = NA", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color="#666666")
    fig.colorbar(image, cax=cax, label="Mean unique-CpG methylation ratio" if raw else "DMR-wise Z-score")
    fig.savefig(output, dpi=dpi); plt.close(fig)


def stage_plots(args, matrix_summary):
    root = args.run_dir / "11_hypo_DMR_heatmaps"; root.mkdir(exist_ok=True)
    summaries = []
    for row in matrix_summary["samples"]:
        sample = row["sample_id"]
        sample_dir = root / "samples" / safe(sample); sample_dir.mkdir(parents=True, exist_ok=True)
        if not row["dmrs"]:
            summaries.append({"sample_id": sample, "status": "no_dmrs", "dmrs": 0, "single_cells": row["single_cells"]}); continue
        cells, matrix_rows, matrix = load_matrix(Path(row["matrix_file"]))
        cell_meta = {r["cell_id"]: r["cell_type"] for r in read_tsv(Path(row["cell_annotations"]))}
        dmr_meta = {r["dmr_id"]: r for r in read_tsv(Path(row["dmr_annotations"]))}
        # Multi-supported merged DMRs are assigned to the supporting hypo type
        # with the lowest observed cell-type mean, matching the reference plot.
        column_types = []
        for j, mr in enumerate(matrix_rows):
            supported = dmr_meta[mr[3]]["hypo_cell_types"].split(",")
            means = []
            for cell_type in supported:
                idx = [i for i, cell in enumerate(cells) if cell_meta[cell] == cell_type]
                value = np.nanmean(matrix[idx, j]) if idx else np.nan
                means.append((not np.isfinite(value), value if np.isfinite(value) else math.inf, cell_type))
            column_types.append(min(means)[2])
        own_dmr_types = set(column_types)
        row_order = sorted(
            (i for i in range(len(cells)) if cell_meta[cells[i]] in own_dmr_types),
            key=lambda i: (cell_meta[cells[i]], cells[i]),
        )
        if not row_order:
            raise ValueError(f"{sample}: no cells belong to a cell type with an assigned hypo-DMR")
        row_types = [cell_meta[cells[i]] for i in row_order]
        col_order = sorted(range(matrix.shape[1]), key=lambda j: (column_types[j], chrom_key(matrix_rows[j][0]), int(matrix_rows[j][1]), int(matrix_rows[j][2])))
        ordered = matrix[np.ix_(row_order, col_order)]
        ordered_col_types = [column_types[j] for j in col_order]
        z, counts, eligible = zscore(ordered, args.zscore_min_observed_cells)
        np.savez_compressed(sample_dir / "heatmap_matrices.npz", mean_ratio=ordered,
                            dmrwise_standard_zscore=z, cells=np.asarray([cells[i] for i in row_order]),
                            dmr_ids=np.asarray([matrix_rows[j][3] for j in col_order]))
        raw_png = sample_dir / "mean_ratio.png"
        z_png = sample_dir / "zscore.png"
        compressed_png = sample_dir / "zscore_colorbar_compressed.png"
        plot_one(ordered, row_types, ordered_col_types, raw_png, f"{sample}: raw mean ratio", 0, 1, True, args.dpi)
        plot_one(np.clip(z, -args.zscore_clip, args.zscore_clip), row_types, ordered_col_types, z_png,
                 f"{sample}: DMR-wise z-score", -args.zscore_clip, args.zscore_clip, False, args.dpi)
        plot_one(z, row_types, ordered_col_types, compressed_png,
                 f"{sample}: DMR-wise z-score (compressed colorbar)", -args.compressed_colorbar_limit,
                 args.compressed_colorbar_limit, False, args.dpi)
        with (sample_dir / "zscore_qc.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["dmr_id", "observed_cells", "eligible"])
            for j, source_j in enumerate(col_order): writer.writerow([matrix_rows[source_j][3], int(counts[j]), "yes" if eligible[j] else "no"])
        summaries.append({"sample_id": sample, "status": "complete", "dmrs": matrix.shape[1],
                          "input_single_cells": matrix.shape[0], "plotted_single_cells": len(row_order),
                          "cell_types_with_own_hypo_dmrs": len(own_dmr_types),
                          "zscore_eligible_dmrs": int(eligible.sum()),
                          "mean_ratio_plot": str(raw_png), "zscore_plot": str(z_png),
                          "zscore_compressed_colorbar_plot": str(compressed_png)})
    summary = {"status": "complete", "zscore_axis": "each DMR across observed single cells",
               "cell_row_rule": "retain only cell types assigned at least one own hypo-DMR",
               "zscore_min_observed_cells": args.zscore_min_observed_cells,
               "ordinary_zscore_display_clip": [-args.zscore_clip, args.zscore_clip],
               "compressed_colorbar_limits": [-args.compressed_colorbar_limit, args.compressed_colorbar_limit],
               "compressed_plot_modifies_zscore_values": False, "samples": summaries}
    (root / "plot_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main():
    args = args_parser(); args.run_dir = args.run_dir.resolve()
    for value in (args.top_per_cell_type, args.min_cells, args.threads, args.matrix_workers, args.zscore_min_observed_cells, args.dpi):
        if value < 1: raise ValueError("Integer parameters must be positive")
    comparisons = read_tsv(args.run_dir / "07_methdiff" / "pairwise_summary.tsv")
    if any(r["status"] == "failed" for r in comparisons):
        raise RuntimeError("Hard pairwise meth-diff failures must be resolved before downstream processing")
    for name in ("08_rawp_fallback", "09_top200_hypo_DMRs", "10_single_cell_DMR_matrix", "11_hypo_DMR_heatmaps"):
        if (args.run_dir / name).exists(): raise FileExistsError(args.run_dir / name)
    fallback_paths, fallback = stage_fallback(args, comparisons)
    selection = stage_select(args, comparisons, fallback_paths)
    matrix = stage_matrix(args, selection)
    plots = stage_plots(args, matrix)
    summary = {"status": "complete", "fallback": fallback, "selection": selection, "matrix": matrix, "plots": plots}
    (args.run_dir / "hypo_dmr_heatmap_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
