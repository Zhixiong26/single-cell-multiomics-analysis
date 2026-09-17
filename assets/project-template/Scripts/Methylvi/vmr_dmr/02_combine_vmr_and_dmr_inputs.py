#!/usr/bin/env python3
"""Append pooled DMR count features that do not overlap a selected VMR input."""
from __future__ import annotations

import argparse
import json
from bisect import bisect_right
from pathlib import Path

import anndata as ad
import mudata
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vmr-input", type=Path, required=True)
    parser.add_argument("--dmr-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    vmr = mudata.read_h5mu(args.vmr_input)["mCG"]
    dmr = mudata.read_h5mu(args.dmr_input)["mCG"]
    if not vmr.obs_names.equals(dmr.obs_names):
        raise ValueError("VMR and DMR cell order differs")
    required = {"chrom", "start", "end"}
    if not required.issubset(vmr.var) or not required.issubset(dmr.var):
        raise KeyError("Both inputs require chrom/start/end feature metadata")
    lookup = {}
    for chrom, frame in vmr.var.groupby("chrom", sort=False):
        intervals = sorted(zip(frame.start.astype(int), frame.end.astype(int)))
        lookup[str(chrom)] = ([x[0] for x in intervals], [x[1] for x in intervals])

    keep = []
    for index, row in enumerate(dmr.var.itertuples()):
        starts_ends = lookup.get(str(row.chrom))
        overlaps = False
        if starts_ends is not None:
            starts, ends = starts_ends
            pos = bisect_right(starts, int(row.end) - 1) - 1
            overlaps = pos >= 0 and ends[pos] > int(row.start)
        if not overlaps:
            keep.append(index)
    dmr_keep = dmr[:, keep].copy()
    if dmr_keep.n_vars == 0:
        raise RuntimeError("No non-overlapping DMR features remain")
    if vmr.var_names.intersection(dmr_keep.var_names).size:
        raise ValueError("VMR and DMR feature IDs collide")
    vmr = vmr.copy()
    vmr.var["feature_source"] = "VMR"
    dmr_keep.var["feature_source"] = "pooled_DMR"
    combined = ad.concat([vmr, dmr_keep], axis=1, join="outer", merge="first")
    combined.var["selection_rank"] = np.arange(1, combined.n_vars + 1, dtype=np.int64)
    for layer in ("mc", "cov"):
        if layer not in combined.layers or not np.issubdtype(combined.layers[layer].dtype, np.integer):
            raise ValueError(f"Combined {layer} layer is absent or non-integer")
    if np.any(combined.layers["mc"] > combined.layers["cov"]):
        raise ValueError("Combined counts contain mc > cov")
    output = mudata.MuData({"mCG": combined})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp.h5mu")
    temporary.unlink(missing_ok=True)
    output.write_h5mu(temporary, compression="gzip")
    temporary.replace(args.output)
    summary = {
        "status": "complete", "vmr_input": str(args.vmr_input.resolve()),
        "dmr_input": str(args.dmr_input.resolve()), "cells": combined.n_obs,
        "vmr_features": vmr.n_vars, "source_merged_dmrs": dmr.n_vars,
        "dmrs_overlapping_vmr_removed": dmr.n_vars - dmr_keep.n_vars,
        "nonoverlapping_dmrs_added": dmr_keep.n_vars,
        "combined_features": combined.n_vars,
        "overlap_rule": "retain VMR; append only DMR intervals with zero base overlap",
    }
    args.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
