#!/usr/bin/env python3
"""Export per-cell sequencing-depth and mCG summaries from a MethylVI H5MU."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-cells", type=int, default=128)
    args = parser.parse_args()
    if not args.input.is_file() or args.chunk_cells < 1:
        raise ValueError("Valid input and positive chunk-cells are required")
    mdata = mudata.read_h5mu(args.input, backed="r")
    adata = mdata["mCG"]
    mc_layer, cov_layer = adata.layers["mc"], adata.layers["cov"]
    depth = np.zeros(adata.n_obs, dtype=np.uint64)
    total_mc = np.zeros(adata.n_obs, dtype=np.uint64)
    mean_mcg = np.full(adata.n_obs, np.nan, dtype=float)
    for start in range(0, adata.n_obs, args.chunk_cells):
        stop = min(start + args.chunk_cells, adata.n_obs)
        mc = np.asarray(mc_layer[start:stop, :])
        cov = np.asarray(cov_layer[start:stop, :])
        if np.any(mc > cov):
            raise ValueError("Invalid input: mc > cov")
        depth[start:stop] = cov.sum(axis=1, dtype=np.uint64)
        total_mc[start:stop] = mc.sum(axis=1, dtype=np.uint64)
        fractions = np.divide(mc, cov, out=np.full(mc.shape, np.nan, dtype=float), where=cov > 0)
        mean_mcg[start:stop] = np.nanmean(fractions, axis=1)
    obs = adata.obs.copy()
    obs["sequencing_depth"] = depth
    obs["overall_mCG"] = np.divide(total_mc, depth, out=np.full(adata.n_obs, np.nan), where=depth > 0)
    obs["mean_mCG"] = mean_mcg
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "cell_methylation_qc.tsv.gz"
    obs.to_csv(table_path, sep="\t", compression="gzip")
    group_key = "sample_id" if "sample_id" in obs else ("condition" if "condition" in obs else None)
    figures = []
    for metric, label in (
        ("sequencing_depth", "Sequencing depth (summed VMR/bin coverage)"),
        ("overall_mCG", "Overall mCG = sum(mc) / sum(cov)"),
        ("mean_mCG", "Mean per-feature mCG"),
    ):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        if group_key:
            groups = [frame[metric].dropna().to_numpy() for _name, frame in obs.groupby(group_key, observed=True)]
            labels = [str(name) for name, _frame in obs.groupby(group_key, observed=True)]
            ax.boxplot(groups, tick_labels=labels, showfliers=False)
            ax.set_xlabel(group_key)
        else:
            ax.hist(obs[metric].dropna().to_numpy(), bins=60)
        ax.set_ylabel(label)
        ax.set_title(metric)
        fig.tight_layout()
        path = args.output_dir / f"{metric}.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        figures.append(str(path))
    summary = {
        "input": str(args.input.resolve()), "cells": int(adata.n_obs),
        "features": int(adata.n_vars), "table": str(table_path), "figures": figures,
    }
    (args.output_dir / "methylation_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
