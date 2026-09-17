#!/usr/bin/env python3
"""Rebuild integer mCG mc/cov for cells and 5-kb bins retained by ALLCools."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gc
import gzip
import hashlib
import json
import os
import re
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd


_LOOKUP: dict[str, dict[int, tuple[int, int, int]]] | None = None
_N_FEATURES = 0
_BIN_SIZE = 5000


def region_table(var: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    lower = {str(column).lower(): str(column) for column in var.columns}
    chrom_column = next((lower[x] for x in ("chrom", "chr", "chromosome", "chrom5k_chrom") if x in lower), None)
    start_column = next((lower[x] for x in ("start", "chrom5k_start") if x in lower), None)
    end_column = next((lower[x] for x in ("end", "chrom5k_end") if x in lower), None)
    if chrom_column and start_column and end_column:
        regions = pd.DataFrame(
            {
                "chrom": var[chrom_column].astype(str).to_numpy(),
                "start": pd.to_numeric(var[start_column], errors="raise").to_numpy(),
                "end": pd.to_numeric(var[end_column], errors="raise").to_numpy(),
            },
            index=var.index.astype(str),
        )
    else:
        patterns = (
            re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)-(?P<end>\d+)$"),
            re.compile(r"^(?P<chrom>chr.+?)[_-](?P<start>\d+)[_-](?P<end>\d+)$"),
        )
        parsed = []
        for feature in var.index.astype(str):
            match = next((pattern.match(feature) for pattern in patterns if pattern.match(feature)), None)
            if match is None:
                raise ValueError(f"Cannot parse ALLCools feature coordinate: {feature!r}")
            parsed.append((match.group("chrom"), int(match.group("start")), int(match.group("end"))))
        regions = pd.DataFrame(parsed, columns=["chrom", "start", "end"], index=var.index.astype(str))
    regions[["start", "end"]] = regions[["start", "end"]].astype(np.int64)
    if (regions["start"] < 0).any() or (regions["end"] <= regions["start"]).any():
        raise ValueError("Invalid retained genomic interval")
    if (regions["start"] % bin_size != 0).any() or ((regions["end"] - regions["start"]) > bin_size).any():
        raise ValueError(f"Retained regions are not aligned to {bin_size:,}-bp bins")
    regions["bin"] = regions["start"] // bin_size
    if regions.duplicated(["chrom", "bin"]).any():
        raise ValueError("Duplicate chromosome/bin coordinates in retained features")
    regions["feature_index"] = np.arange(len(regions), dtype=np.int64)
    return regions


def make_lookup(regions: pd.DataFrame) -> dict[str, dict[int, tuple[int, int, int]]]:
    lookup: dict[str, dict[int, tuple[int, int, int]]] = {}
    for row in regions.itertuples(index=False):
        lookup.setdefault(str(row.chrom), {})[int(row.bin)] = (
            int(row.feature_index), int(row.start), int(row.end)
        )
    return lookup


def init_worker(lookup: dict[str, dict[int, tuple[int, int, int]]], n_features: int, bin_size: int) -> None:
    global _LOOKUP, _N_FEATURES, _BIN_SIZE
    _LOOKUP, _N_FEATURES, _BIN_SIZE = lookup, n_features, bin_size


def checkpoint_valid(path: Path, expected_cell: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with np.load(path, allow_pickle=False) as row:
            return (
                str(row["cell_id"].item()) == expected_cell
                and row["indices"].ndim == row["mc"].ndim == row["cov"].ndim == 1
                and len(row["indices"]) == len(row["mc"]) == len(row["cov"])
                and np.all(row["mc"] <= row["cov"])
            )
    except (OSError, ValueError, KeyError):
        return False


def build_row(task: tuple[int, str, str, str]) -> dict[str, object]:
    row_index, cell, allc_string, checkpoint_string = task
    checkpoint = Path(checkpoint_string)
    if checkpoint_valid(checkpoint, cell):
        return {"row": row_index, "status": "reused"}
    if _LOOKUP is None:
        raise RuntimeError("Worker lookup was not initialized")
    mc = np.zeros(_N_FEATURES, dtype=np.uint64)
    cov = np.zeros(_N_FEATURES, dtype=np.uint64)
    with gzip.open(allc_string, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{allc_string}:{line_number}: expected >=6 columns")
            if not fields[3].upper().startswith("CG"):
                continue
            try:
                position0, methylated, coverage = int(fields[1]) - 1, int(fields[4]), int(fields[5])
            except ValueError as exc:
                raise ValueError(f"{allc_string}:{line_number}: non-integer position/count") from exc
            if position0 < 0 or methylated < 0 or coverage < methylated:
                raise ValueError(f"{allc_string}:{line_number}: invalid mc/cov")
            chrom_bins = _LOOKUP.get(fields[0])
            if chrom_bins is None:
                continue
            region = chrom_bins.get(position0 // _BIN_SIZE)
            if region is None:
                continue
            feature_index, start, end = region
            if start <= position0 < end:
                mc[feature_index] += methylated
                cov[feature_index] += coverage
    nonzero = np.flatnonzero(cov)
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temporary = checkpoint.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary, cell_id=np.asarray(cell), indices=nonzero.astype(np.int32),
        mc=mc[nonzero], cov=cov[nonzero],
        max_mc=np.asarray(mc.max(initial=0), dtype=np.uint64),
        max_cov=np.asarray(cov.max(initial=0), dtype=np.uint64),
    )
    temporary.replace(checkpoint)
    return {"row": row_index, "status": "built"}


def digest_values(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode())
        digest.update(b"\0")
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=Path(os.environ["SCMO_ALLCOOLS_H5AD"]))
    parser.add_argument("--allc-table", type=Path, default=Path(os.environ["SCMO_ALLC_TABLE"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["SCMO_MVI_INPUT"]))
    parser.add_argument("--work-dir", type=Path, default=Path(os.environ["SCMO_MVI_ROOT"]) / "count_rows")
    parser.add_argument("--threads", type=int, default=int(os.environ["SCMO_THREADS"]))
    parser.add_argument("--bin-size", type=int, default=int(os.environ["SCMO_BIN_SIZE"]))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.threads < 1 or args.bin_size < 1:
        raise ValueError("threads and bin-size must be positive")
    source = ad.read_h5ad(args.h5ad, backed="r")
    cells = source.obs_names.astype(str).tolist()
    regions = region_table(source.var, args.bin_size)
    features = regions.index.astype(str).tolist()
    feature_coordinates = [
        f"{feature}\t{row.chrom}\t{int(row.start)}\t{int(row.end)}"
        for feature, row in regions.iterrows()
    ]
    manifest = {
        "h5ad": str(args.h5ad.resolve()), "cells": len(cells), "features": len(features),
        "cell_sha256": digest_values(cells),
        "feature_sha256": digest_values(feature_coordinates),
        "obs_sha256": hashlib.sha256(source.obs.to_csv(sep="\t").encode()).hexdigest(),
        "bin_size": args.bin_size, "context": "CGN",
    }
    args.work_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.work_dir / "manifest.json"
    rows_present = list(args.work_dir.glob("*.npz"))
    if manifest_path.exists():
        if json.loads(manifest_path.read_text()) != manifest:
            raise RuntimeError(f"Checkpoint manifest differs: use a new work directory: {args.work_dir}")
    elif rows_present:
        raise RuntimeError(f"Unversioned checkpoints found in {args.work_dir}")
    else:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    shape = (len(cells), len(regions))
    if args.output.is_file() and args.output.stat().st_size > 0:
        existing = mudata.read_h5mu(args.output, backed="r")
        compatible = (
            "mCG" in existing.mod
            and existing["mCG"].shape == shape
            and {"mc", "cov"}.issubset(existing["mCG"].layers)
            and digest_values(existing["mCG"].obs_names.astype(str).tolist()) == manifest["cell_sha256"]
            and digest_values([
                f"{feature}\t{row.chrom}\t{int(row.start)}\t{int(row.end)}"
                for feature, row in existing["mCG"].var.iterrows()
            ]) == manifest["feature_sha256"]
        )
        if getattr(existing, "file", None) is not None:
            existing.file.close()
        if not compatible:
            raise RuntimeError(f"Existing H5MU is incompatible with the current manifest: {args.output}")
        source.file.close()
        print(f"Existing compatible H5MU detected: {args.output}")
        return

    table = pd.read_csv(args.allc_table, sep="\t", header=None, names=["cell_id", "path"], dtype=str)
    if table["cell_id"].duplicated().any():
        raise ValueError("ALLC table contains duplicate cell IDs")
    allc_index = table.set_index("cell_id")["path"].to_dict()
    missing = [cell for cell in cells if cell not in allc_index]
    if missing:
        raise FileNotFoundError(f"No ALLC matched {len(missing)} selected cells; first: {missing[:5]}")
    row_paths, tasks = [], []
    for row_index, cell in enumerate(cells):
        token = hashlib.sha1(cell.encode(), usedforsecurity=False).hexdigest()[:10]
        row_path = args.work_dir / f"{row_index:06d}.{token}.npz"
        row_paths.append(row_path)
        tasks.append((row_index, cell, allc_index[cell], str(row_path)))

    lookup = make_lookup(regions)
    built = reused = 0
    with cf.ProcessPoolExecutor(
        max_workers=args.threads, initializer=init_worker,
        initargs=(lookup, len(regions), args.bin_size),
    ) as executor:
        for completed, result in enumerate(executor.map(build_row, tasks), start=1):
            built += result["status"] == "built"
            reused += result["status"] == "reused"
            if completed % 25 == 0 or completed == len(tasks):
                print(f"count rows {completed:,}/{len(tasks):,} (built={built:,}, reused={reused:,})", flush=True)

    max_cov = 0
    for cell, path in zip(cells, row_paths):
        if not checkpoint_valid(path, cell):
            raise RuntimeError(f"Invalid count checkpoint: {path}")
        with np.load(path, allow_pickle=False) as row:
            max_cov = max(max_cov, int(row["max_cov"]))
    dtype = np.dtype("uint16" if max_cov <= np.iinfo(np.uint16).max else "uint32")
    mmap_dir = args.work_dir / "assembly"
    mmap_dir.mkdir(exist_ok=True)
    mc = np.memmap(mmap_dir / f"mc.{dtype.name}.mmap", mode="w+", dtype=dtype, shape=shape)
    cov = np.memmap(mmap_dir / f"cov.{dtype.name}.mmap", mode="w+", dtype=dtype, shape=shape)
    mc[:] = 0
    cov[:] = 0
    for row_index, path in enumerate(row_paths):
        with np.load(path, allow_pickle=False) as row:
            indices = row["indices"]
            mc[row_index, indices] = row["mc"].astype(dtype, copy=False)
            cov[row_index, indices] = row["cov"].astype(dtype, copy=False)
        if (row_index + 1) % 100 == 0 or row_index + 1 == len(row_paths):
            print(f"assembled {row_index + 1:,}/{len(row_paths):,}", flush=True)
    mc.flush()
    cov.flush()

    obs = source.obs.copy()
    var = source.var.copy()
    var["chrom"] = regions["chrom"].to_numpy()
    var["start"] = regions["start"].to_numpy()
    var["end"] = regions["end"].to_numpy()
    counts = ad.AnnData(X=None, obs=obs, var=var)
    counts.layers["mc"] = np.asarray(mc)
    counts.layers["cov"] = np.asarray(cov)
    output = mudata.MuData({"mCG": counts})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp.h5mu")
    temporary.unlink(missing_ok=True)
    output.write_h5mu(temporary, compression="gzip")
    temporary.replace(args.output)
    source.file.close()

    reopened = mudata.read_h5mu(args.output, backed="r")
    if reopened["mCG"].shape != shape or not {"mc", "cov"}.issubset(reopened["mCG"].layers):
        raise RuntimeError("Written H5MU failed shape/layer verification")
    if getattr(reopened, "file", None) is not None:
        reopened.file.close()
    summary = {
        **manifest, "output": str(args.output), "dtype": dtype.name,
        "maximum_cov_per_cell_bin": max_cov, "count_rows_built": built,
        "count_rows_reused": reused,
    }
    (args.output.parent / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    del output, counts, mc, cov
    gc.collect()
    for path in mmap_dir.glob("*.mmap"):
        path.unlink()
    if not any(mmap_dir.iterdir()):
        mmap_dir.rmdir()


if __name__ == "__main__":
    main()
