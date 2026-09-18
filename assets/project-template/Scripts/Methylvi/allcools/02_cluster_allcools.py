#!/usr/bin/env python3
"""Apply the yuanpei ALLCools mCG 5-kb clustering logic."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from ALLCools.clustering import (
    ConsensusClustering,
    binarize_matrix,
    lsi,
    significant_pc_test,
    tsne,
)
from ALLCools.mcds import MCDS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcds", type=Path, default=Path(os.environ["SCMO_MCDS"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["SCMO_ALLCOOLS_H5AD"]))
    parser.add_argument("--annotation", type=Path, default=Path(os.environ["SCMO_ANNOTATION"]))
    parser.add_argument("--blacklist", type=Path, default=Path(os.environ["SCMO_BLACKLIST"]))
    parser.add_argument("--blacklist-md5", default=os.environ["SCMO_BLACKLIST_MD5"])
    parser.add_argument(
        "--blacklist-fraction", type=float,
        default=float(os.environ["SCMO_BLACKLIST_FRACTION"]),
    )
    parser.add_argument(
        "--target-features", type=int, default=int(os.environ["SCMO_TARGET_FEATURES"]),
        help="Select exactly this many top 5-kb bins by current-cell hypo prevalence",
    )
    parser.add_argument("--threads", type=int, default=int(os.environ["SCMO_THREADS"]))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not 0 < args.blacklist_fraction <= 1:
        raise ValueError("blacklist fraction must be in (0, 1]")
    if args.target_features < 2:
        raise ValueError("target-features must be at least 2")
    if not args.blacklist.is_file():
        raise FileNotFoundError(args.blacklist)
    blacklist_md5 = hashlib.md5(args.blacklist.read_bytes()).hexdigest()
    if args.blacklist_md5 and blacklist_md5.lower() != args.blacklist_md5.lower():
        raise ValueError(
            f"Blacklist MD5 mismatch: expected={args.blacklist_md5} observed={blacklist_md5}"
        )
    warnings.filterwarnings("ignore", category=FutureWarning)
    mcds = MCDS.open(str(args.mcds), var_dim="chrom5k")
    bins_before_blacklist = int(mcds.get_index("chrom5k").size)
    mcds = mcds.remove_black_list_region(
        black_list_path=str(args.blacklist), f=args.blacklist_fraction,
    )
    bins_after_blacklist = int(mcds.get_index("chrom5k").size)
    print(
        f"Blacklist removed {bins_before_blacklist - bins_after_blacklist:,} / "
        f"{bins_before_blacklist:,} 5-kb bins (overlap >= {args.blacklist_fraction:g})",
        flush=True,
    )
    adata = mcds.get_score_adata(mc_type=os.environ["SCMO_MC_CONTEXT"], quant_type="hypo-score")
    initial_shape = [int(adata.n_obs), int(adata.n_vars)]
    print(f"Initial matrix: {adata.n_obs:,} cells x {adata.n_vars:,} bins", flush=True)

    binarize_matrix(adata, cutoff=float(os.environ["SCMO_BINARIZE_CUTOFF"]))
    prevalence_counts = np.asarray(adata.X.sum(axis=0)).ravel().astype(np.int64)
    feature_names = adata.var_names.astype(str).to_numpy()
    eligible = np.flatnonzero(prevalence_counts > 0)
    if len(eligible) < args.target_features:
        raise RuntimeError(
            f"Only {len(eligible):,} bins have nonzero hypo prevalence; "
            f"cannot select {args.target_features:,}"
        )
    ranked = eligible[np.lexsort((feature_names[eligible], -prevalence_counts[eligible]))]
    selected = ranked[: args.target_features]
    selected_counts = prevalence_counts[selected]
    boundary_count = int(selected_counts[-1])
    boundary_ties = int((prevalence_counts == boundary_count).sum())
    adata = adata[:, selected].copy()
    adata.var["hypo_cells"] = selected_counts
    adata.var["hypo_percent"] = selected_counts / adata.n_obs * 100.0
    adata.var["selection_rank"] = np.arange(1, adata.n_vars + 1, dtype=np.int64)
    effective_hypo_percent = boundary_count / adata.n_obs * 100.0
    feature_filter_summary = {
        "cells": int(adata.n_obs),
        "bins_after_blacklist": bins_after_blacklist,
        "nonzero_hypo_bins": int(len(eligible)),
        "target_features": args.target_features,
        "retained_features": int(adata.n_vars),
        "boundary_hypo_cells": boundary_count,
        "effective_hypo_percent": effective_hypo_percent,
        "MVI_HYPO_PERCENT": effective_hypo_percent,
        "boundary_tied_bins": boundary_ties,
        "ranking": "hypo_cells descending, feature_id ascending",
    }
    if adata.n_vars != args.target_features:
        raise RuntimeError("Feature target hard check failed")
    (args.output.parent / "feature_filter_summary.json").write_text(
        json.dumps(feature_filter_summary, indent=2) + "\n"
    )
    filtered_shape = [int(adata.n_obs), int(adata.n_vars)]
    print(f"Filtered matrix: {adata.n_obs:,} cells x {adata.n_vars:,} bins", flush=True)
    seed = int(os.environ["SCMO_SEED"])
    requested_components = int(os.environ["SCMO_LSI_COMPONENTS"])
    # scipy.sparse.linalg.svds (used by ARPACK) requires
    # 0 < k < min(matrix.shape). Keep the full-run value at 100 while making
    # small bounded validation datasets valid.
    lsi_components = min(requested_components, adata.n_obs - 1, adata.n_vars - 1)
    if lsi_components < 1:
        raise ValueError(
            f"LSI requires at least 2 cells and 2 retained bins; got {adata.shape}"
        )
    print(
        f"LSI components: {lsi_components} (requested {requested_components})",
        flush=True,
    )
    lsi(
        adata, n_components=lsi_components, algorithm="arpack",
        obsm="X_pca", random_state=seed,
    )
    n_components = significant_pc_test(
        adata, p_cutoff=float(os.environ["SCMO_LSI_P_CUTOFF"]), update=True
    )
    print(f"Significant LSI components: {n_components}", flush=True)

    neighbors = int(os.environ["SCMO_ALLCOOLS_NEIGHBORS"])
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=neighbors, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=float(os.environ["SCMO_ALLCOOLS_LEIDEN_RESOLUTION"]),
        random_state=seed,
        key_added="leiden",
    )
    tsne(
        adata, obsm="X_pca", metric="euclidean", exaggeration=-1,
        perplexity=30, n_jobs=args.threads,
    )
    sc.tl.umap(adata, random_state=seed)
    consensus = ConsensusClustering(
        model=None, n_neighbors=neighbors, metric="euclidean", min_cluster_size=10,
        leiden_repeats=int(os.environ["SCMO_CONSENSUS_LEIDEN_REPEATS"]),
        leiden_resolution=float(os.environ["SCMO_CONSENSUS_LEIDEN_RESOLUTION"]),
        consensus_rate=0.5, random_state=seed, train_frac=0.5, train_max_n=500,
        max_iter=20, n_jobs=args.threads,
    )
    consensus.fit_predict(adata.obsm["X_pca"])
    adata.obs["L1"] = pd.Categorical(np.asarray(consensus.label).astype(str))
    adata.obs["L1_proba"] = np.asarray(consensus.label_proba, dtype=float)

    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    if annotation["cell_id"].duplicated().any():
        raise ValueError("Annotation contains duplicate cell_id values")
    annotation = annotation.set_index("cell_id")
    aligned = annotation.reindex(adata.obs_names)
    celltype_column = "cell_type" if "cell_type" in aligned.columns else "manual_celltype"
    adata.obs["cell_type"] = aligned[celltype_column].fillna("Unknown").to_numpy()
    adata.obs["manual_celltype"] = adata.obs["cell_type"].to_numpy()
    adata.obs["sample_id"] = (
        aligned["sample"].fillna(adata.obs_names.to_series().str.split("_", n=1).str[0]).to_numpy()
        if "sample" in aligned.columns else adata.obs_names.to_series().str.split("_", n=1).str[0].to_numpy()
    )
    adata.obs["condition"] = (
        aligned["cohort"].fillna(adata.obs["sample_id"]).to_numpy()
        if "cohort" in aligned.columns else adata.obs["sample_id"].to_numpy()
    )
    adata.obs["cohort"] = adata.obs["condition"].to_numpy()
    adata.write_h5ad(args.output, compression="gzip")
    adata.obs.to_csv(args.output.parent / "cell_clusters.csv.gz")
    for basis in ("tsne", "umap"):
        sc.pl.embedding(
            adata, basis=basis,
            color=["L1", "L1_proba", "cohort", "manual_celltype"],
            show=False, wspace=0.35,
        )
        plt.savefig(args.output.parent / f"allcools_{basis}.png", dpi=300, bbox_inches="tight")
        plt.close("all")
    summary = {
        "initial_shape": initial_shape, "filtered_shape": filtered_shape,
        "blacklist": str(args.blacklist.resolve()), "blacklist_md5": blacklist_md5,
        "blacklist_fraction": args.blacklist_fraction,
        "bins_before_blacklist": bins_before_blacklist,
        "bins_after_blacklist": bins_after_blacklist,
        "blacklist_removed_bins": bins_before_blacklist - bins_after_blacklist,
        "target_features": args.target_features,
        "effective_hypo_percent": effective_hypo_percent,
        "requested_lsi_components": requested_components,
        "computed_lsi_components": lsi_components,
        "significant_lsi_components": int(n_components), "neighbors": neighbors, "seed": seed,
    }
    (args.output.parent / "cluster_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
