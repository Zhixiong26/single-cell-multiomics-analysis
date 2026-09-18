#!/usr/bin/env python3
"""Aggregate integer mCG mc/cov counts over non-overlapping VMRs for MethylVI."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gc
import gzip
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd


_REGIONS = None
_N_REGIONS = 0


def file_sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_regions(path: Path) -> tuple[pd.DataFrame, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    table = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "peak_var", "n_cpg", "n_obs_cells", "vmr_id"],
        dtype={"chrom": str, "vmr_id": str},
    )
    if table.empty or table["vmr_id"].duplicated().any():
        raise ValueError("VMR BED must contain non-empty, unique region IDs")
    table[["start", "end"]] = table[["start", "end"]].astype(np.int64)
    table[["peak_var", "n_cpg", "n_obs_cells"]] = table[["peak_var", "n_cpg", "n_obs_cells"]].apply(pd.to_numeric, errors="raise")
    if (table.start < 0).any() or (table.end <= table.start).any():
        raise ValueError("Invalid VMR coordinates")
    table = table.set_index("vmr_id", drop=True)
    lookup = {}
    for chrom, frame in table.groupby("chrom", sort=False):
        indices = table.index.get_indexer(frame.index).astype(np.int32)
        starts = frame.start.to_numpy(np.int64)
        ends = frame.end.to_numpy(np.int64)
        if np.any(starts[1:] < ends[:-1]):
            raise ValueError(f"Overlapping VMRs on {chrom} are not supported")
        lookup[str(chrom)] = (starts, ends, indices)
    return table, lookup


def init_worker(regions, n_regions: int) -> None:
    global _REGIONS, _N_REGIONS
    _REGIONS, _N_REGIONS = regions, n_regions


def checkpoint_valid(path: Path, cell: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with np.load(path, allow_pickle=False) as row:
            indices, mc, cov = row["indices"], row["mc"], row["cov"]
            return (
                str(row["cell_id"].item()) == cell and indices.ndim == mc.ndim == cov.ndim == 1
                and len(indices) == len(mc) == len(cov) and np.all(indices < _N_REGIONS)
                and np.all(mc <= cov)
            )
    except (OSError, ValueError, KeyError):
        return False


def build_row(task: tuple[int, str, str, str]) -> dict[str, object]:
    row_index, cell, allc_string, checkpoint_string = task
    checkpoint = Path(checkpoint_string)
    if checkpoint_valid(checkpoint, cell):
        return {"row": row_index, "status": "reused"}
    if _REGIONS is None:
        raise RuntimeError("Worker region lookup was not initialized")
    mc = np.zeros(_N_REGIONS, dtype=np.uint64)
    cov = np.zeros(_N_REGIONS, dtype=np.uint64)
    pointers: dict[str, int] = {}
    with gzip.open(allc_string, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip().split("\t")
            if len(fields) < 6:
                raise ValueError(f"{allc_string}:{line_number}: expected >=6 columns")
            if not fields[3].upper().startswith("CG"):
                continue
            region = _REGIONS.get(fields[0])
            if region is None:
                continue
            try:
                position0, methylated, coverage = int(fields[1]) - 1, int(fields[4]), int(fields[5])
            except ValueError as exc:
                raise ValueError(f"{allc_string}:{line_number}: non-integer position/count") from exc
            if position0 < 0 or methylated < 0 or coverage < methylated:
                raise ValueError(f"{allc_string}:{line_number}: invalid mc/cov")
            starts, ends, indices = region
            pointer = pointers.get(fields[0], 0)
            while pointer < len(starts) and position0 >= ends[pointer]:
                pointer += 1
            pointers[fields[0]] = pointer
            if pointer < len(starts) and starts[pointer] <= position0 < ends[pointer]:
                index = indices[pointer]
                mc[index] += methylated
                cov[index] += coverage
    nonzero = np.flatnonzero(cov)
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temporary = checkpoint.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary, cell_id=np.asarray(cell), indices=nonzero.astype(np.int32),
        mc=mc[nonzero], cov=cov[nonzero], max_cov=np.asarray(cov.max(initial=0), dtype=np.uint64),
    )
    temporary.replace(checkpoint)
    return {"row": row_index, "status": "built"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, default=Path(os.environ["VMR_FILTERED_BED"]))
    parser.add_argument("--allc-table", type=Path, default=Path(os.environ["VMR_ALLC_TABLE"]))
    parser.add_argument("--annotation", type=Path, default=Path(os.environ["VMR_ANNOTATION"]))
    parser.add_argument("--work-dir", type=Path, default=Path(os.environ["VMR_COUNT_ROWS"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["VMR_MVI_INPUT"]))
    parser.add_argument("--threads", type=int, default=int(os.environ["VMR_THREADS"]))
    parser.add_argument("--min-covered-percent", type=float, default=float(os.environ["VMR_MIN_COVERED_PERCENT"]))
    parser.add_argument("--target-features", type=int, default=int(os.environ["VMR_TARGET_FEATURES"]))
    parser.add_argument("--max-cells", type=int, default=0, help="Use only the first N cells for a smoke test")
    parser.add_argument("--mc-context", default=os.environ.get("VMR_MC_CONTEXT", os.environ.get("SCMO_MC_CONTEXT", "CGN")))
    args = parser.parse_args()
    if args.threads < 1 or not 0 <= args.min_covered_percent <= 100 or args.max_cells < 0 or args.target_features < 2:
        raise ValueError("threads must be positive; min-covered-percent must be in [0,100]; max-cells non-negative")
    for path in (args.bed, args.allc_table, args.annotation):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    regions, lookup = load_regions(args.bed)
    if len(regions) < args.target_features:
        raise RuntimeError(f"Only {len(regions):,} source VMRs; cannot select {args.target_features:,}")
    # Initialize the parent process as well as ProcessPool workers so checkpoint
    # validation uses the same feature bound during assembly.
    init_worker(lookup, len(regions))
    allc_table = pd.read_csv(args.allc_table, sep="\t", header=None, names=["cell_id", "path"], dtype=str)
    if allc_table.empty or allc_table.cell_id.duplicated().any():
        raise ValueError("ALLC table is empty or contains duplicate cell IDs")
    if args.max_cells:
        allc_table = allc_table.iloc[: args.max_cells].copy()
    cells = allc_table.cell_id.tolist()
    configured_floor_cells = int(np.ceil(len(cells) * args.min_covered_percent / 100.0))
    args.work_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "bed": str(args.bed.resolve()), "bed_sha256": file_sha256(args.bed),
        "allc_table_sha256": file_sha256(args.allc_table), "annotation_sha256": file_sha256(args.annotation),
        "cells": len(cells), "source_vmrs": len(regions), "context": args.mc_context,
    }
    manifest_path = args.work_dir / "manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise RuntimeError(f"Checkpoint manifest differs; use a new work directory: {args.work_dir}")
    if not manifest_path.exists():
        if any(args.work_dir.glob("*.npz")):
            raise RuntimeError(f"Unversioned checkpoints found: {args.work_dir}")
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    tasks, row_paths = [], []
    for row_index, row in enumerate(allc_table.itertuples(index=False)):
        allc = Path(row.path)
        if not allc.is_file() or allc.stat().st_size == 0:
            raise FileNotFoundError(allc)
        token = hashlib.sha1(row.cell_id.encode(), usedforsecurity=False).hexdigest()[:10]
        checkpoint = args.work_dir / f"{row_index:06d}.{token}.npz"
        row_paths.append(checkpoint)
        tasks.append((row_index, row.cell_id, str(allc), str(checkpoint)))
    built = reused = 0
    with cf.ProcessPoolExecutor(
        max_workers=args.threads, initializer=init_worker, initargs=(lookup, len(regions)),
    ) as executor:
        for completed, result in enumerate(executor.map(build_row, tasks), start=1):
            built += result["status"] == "built"
            reused += result["status"] == "reused"
            if completed % 25 == 0 or completed == len(tasks):
                print(f"VMR count rows {completed:,}/{len(tasks):,} (built={built:,}, reused={reused:,})", flush=True)

    covered_cells = np.zeros(len(regions), dtype=np.uint32)
    max_cov = 0
    for cell, checkpoint in zip(cells, row_paths):
        if not checkpoint_valid(checkpoint, cell):
            raise RuntimeError(f"Invalid checkpoint: {checkpoint}")
        with np.load(checkpoint, allow_pickle=False) as row:
            covered_cells[row["indices"]] += 1
            max_cov = max(max_cov, int(row["max_cov"]))
    target_coverage_cutoff = int(np.partition(covered_cells, len(covered_cells) - args.target_features)[len(covered_cells) - args.target_features])
    min_covered_cells = max(configured_floor_cells, target_coverage_cutoff)
    eligible = np.flatnonzero(covered_cells >= min_covered_cells)
    if len(eligible) < args.target_features:
        raise RuntimeError(
            f"Only {len(eligible):,} VMRs pass coverage filtering; cannot select {args.target_features:,}"
        )
    peak = regions["peak_var"].to_numpy(float)
    observed = regions["n_obs_cells"].to_numpy(np.int64)
    cpgs = regions["n_cpg"].to_numpy(np.int64)
    ids = regions.index.astype(str).to_numpy()
    ranked = eligible[np.lexsort((ids[eligible], -cpgs[eligible], -observed[eligible], -peak[eligible]))]
    selected = ranked[: args.target_features]
    retained = len(selected)
    print(
        f"Selected top {retained:,} from {len(eligible):,}/{len(regions):,} coverage-eligible VMRs ",
        f"(>={min_covered_cells} cells; {min_covered_cells / len(cells) * 100.0:.4g}%)",
        flush=True,
    )
    old_to_new = np.full(len(regions), -1, dtype=np.int64)
    old_to_new[selected] = np.arange(retained)
    shape = (len(cells), retained)
    dtype = np.dtype("uint16" if max_cov <= np.iinfo(np.uint16).max else "uint32")
    assembly = args.work_dir / "assembly"
    assembly.mkdir(exist_ok=True)
    mc = np.memmap(assembly / f"mc.{dtype.name}.mmap", mode="w+", dtype=dtype, shape=shape)
    cov = np.memmap(assembly / f"cov.{dtype.name}.mmap", mode="w+", dtype=dtype, shape=shape)
    mc[:] = 0
    cov[:] = 0
    for row_index, checkpoint in enumerate(row_paths):
        with np.load(checkpoint, allow_pickle=False) as row:
            target = old_to_new[row["indices"]]
            selected_mask = target >= 0
            mc[row_index, target[selected_mask]] = row["mc"][selected_mask].astype(dtype, copy=False)
            cov[row_index, target[selected_mask]] = row["cov"][selected_mask].astype(dtype, copy=False)
        if (row_index + 1) % 250 == 0 or row_index + 1 == len(row_paths):
            print(f"Assembled {row_index + 1:,}/{len(row_paths):,} cells", flush=True)
    mc.flush()
    cov.flush()

    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    if "cell_id" not in annotation.columns or not ({"cell_type", "manual_celltype"} & set(annotation.columns)) or annotation.cell_id.duplicated().any():
        raise ValueError("Annotation requires unique cell_id and cell_type/manual_celltype")
    annotation = annotation.set_index("cell_id")
    aligned = annotation.reindex(cells)
    celltype_column = "cell_type" if "cell_type" in aligned.columns else "manual_celltype"
    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    obs["cell_type"] = aligned[celltype_column].fillna("Unknown").to_numpy()
    obs["manual_celltype"] = obs["cell_type"].to_numpy()
    fallback = obs.index.to_series().str.split("_", n=1).str[0]
    obs["sample_id"] = aligned["sample"].fillna(fallback).to_numpy() if "sample" in aligned.columns else fallback.to_numpy()
    obs["condition"] = aligned["cohort"].fillna(obs["sample_id"]).to_numpy() if "cohort" in aligned.columns else obs["sample_id"].to_numpy()
    obs["cohort"] = obs["condition"].to_numpy()
    var = regions.iloc[selected].copy()
    var["length"] = var.end - var.start
    var["covered_cells"] = covered_cells[selected]
    var["selection_rank"] = np.arange(1, retained + 1, dtype=np.int64)
    counts = ad.AnnData(X=None, obs=obs, var=var)
    counts.layers["mc"] = np.asarray(mc)
    counts.layers["cov"] = np.asarray(cov)
    output = mudata.MuData({"mCG": counts})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = args.output.with_suffix(".tmp.h5mu")
    temporary_output.unlink(missing_ok=True)
    output.write_h5mu(temporary_output, compression="gzip")
    temporary_output.replace(args.output)
    summary = {
        **manifest, "output": str(args.output.resolve()),
        "configured_min_covered_percent": args.min_covered_percent,
        "configured_floor_cells": configured_floor_cells,
        "min_covered_cells": min_covered_cells,
        "effective_covered_percent": min_covered_cells / len(cells) * 100.0,
        "coverage_eligible_vmrs": int(len(eligible)),
        "target_features": args.target_features,
        "retained_vmrs": retained,
        "ranking": "peak_var desc, n_obs_cells desc, n_cpg desc, vmr_id asc",
        "dtype": dtype.name, "maximum_cov_per_cell_vmr": max_cov,
        "count_rows_built": built, "count_rows_reused": reused,
        "covered_cells_quantiles": {
            str(q): float(np.quantile(covered_cells, q)) for q in (0, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
        },
    }
    (args.output.parent / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    del output, counts, mc, cov
    gc.collect()
    for mmap in assembly.glob("*.mmap"):
        mmap.unlink()
    if not any(assembly.iterdir()):
        assembly.rmdir()
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
