#!/usr/bin/env python3
"""Create target_weight-supervised UMAPs from the finished MethylVI latent space."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import umap


def parse_weights(raw: str) -> list[float]:
    weights = [float(value) for value in raw.replace(",", " ").split()]
    if not weights or len(weights) != len(set(weights)) or any(not 0 <= x <= 1 for x in weights):
        raise ValueError("weights must be unique values in [0, 1]")
    return weights


def tag(weight: float) -> str:
    return f"{weight:g}".replace(".", "p")


def target_codes(obs: pd.DataFrame, target_key: str) -> tuple[np.ndarray, dict[str, int]]:
    if target_key not in obs:
        raise KeyError(f"Embedding lacks supervised target column: {target_key}")
    labels = obs[target_key].astype("string").fillna("Unknown").str.strip()
    unknown = labels.str.lower().isin({"", "unknown", "unannotated", "na", "nan"}).to_numpy()
    categories = sorted(labels.loc[~unknown].unique())
    if len(categories) < 2:
        raise ValueError("At least two known cell-type labels are required")
    mapping = {label: index for index, label in enumerate(categories)}
    codes = np.full(len(labels), -1, dtype=np.int32)
    for label, code in mapping.items():
        codes[(labels == label).to_numpy()] = code
    return codes, mapping


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(os.environ["SCMO_MVI_RESULTS"]) / "methylvi_embedding.h5ad")
    parser.add_argument("--output", type=Path, default=Path(os.environ["SCMO_MVI_RESULTS"]) / "supervised_umap")
    parser.add_argument("--target-key", default=os.environ["SCMO_SUPERVISED_TARGET_KEY"])
    parser.add_argument("--weights", nargs="+", type=float, default=parse_weights(os.environ["SCMO_SUPERVISED_TARGET_WEIGHTS"]))
    parser.add_argument("--neighbors", type=int, default=int(os.environ["SCMO_SUPERVISED_NEIGHBORS"]))
    parser.add_argument("--min-dist", type=float, default=float(os.environ["SCMO_SUPERVISED_MIN_DIST"]))
    parser.add_argument("--seed", type=int, default=int(os.environ["SCMO_SEED"]))
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    if not args.weights or len(args.weights) != len(set(args.weights)) or any(not 0 <= x <= 1 for x in args.weights):
        raise ValueError("weights must be unique values in [0, 1]")

    embedding = ad.read_h5ad(args.input)
    if "X_methylVI" not in embedding.obsm:
        raise KeyError("Embedding lacks obsm['X_methylVI']")
    latent = np.asarray(embedding.obsm["X_methylVI"], dtype=np.float32)
    if latent.ndim != 2 or latent.shape[0] != embedding.n_obs or not np.isfinite(latent).all():
        raise ValueError("Invalid X_methylVI latent representation")
    labels, label_map = target_codes(embedding.obs, args.target_key)
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    color_keys = [key for key in (args.target_key, os.environ["SCMO_BATCH_KEY"], "L1", "methylVI_leiden") if key in embedding.obs]
    figure_files: list[str] = []
    coordinate_files: dict[str, str] = {}

    for weight in args.weights:
        weight_tag = tag(weight)
        reducer = umap.UMAP(
            n_neighbors=args.neighbors, n_components=2, metric="euclidean", min_dist=args.min_dist,
            target_metric="categorical", target_weight=weight, random_state=args.seed,
            transform_seed=args.seed, n_jobs=1,
        )
        coordinates = reducer.fit_transform(latent, y=labels)
        if coordinates.shape != (embedding.n_obs, 2) or not np.isfinite(coordinates).all():
            raise RuntimeError(f"Invalid coordinates for target_weight={weight:g}")
        basis = f"umap_target_weight_{weight_tag}"
        embedding.obsm[f"X_{basis}"] = np.asarray(coordinates, dtype=np.float32)
        table = embedding.obs.copy()
        table.insert(0, "UMAP1", coordinates[:, 0])
        table.insert(1, "UMAP2", coordinates[:, 1])
        table["target_weight"] = weight
        coordinate_path = output / f"target_weight_{weight_tag}_coordinates.tsv.gz"
        table.to_csv(coordinate_path, sep="\t")
        coordinate_files[f"{weight:g}"] = str(coordinate_path)
        figure_dir = output / f"target_weight_{weight_tag}"
        figure_dir.mkdir(exist_ok=True)
        for color_key in color_keys:
            figure = sc.pl.embedding(
                embedding, basis=basis, color=color_key, show=False, return_fig=True,
                title=f"MethylVI supervised UMAP — {color_key} (target_weight={weight:g})",
            )
            figure_path = figure_dir / f"methylvi_supervised_umap_{color_key}.png"
            figure.savefig(figure_path, dpi=200, bbox_inches="tight")
            figure_files.append(str(figure_path))

    embedding.uns["supervised_umap"] = {
        "source_representation": "X_methylVI", "target_key": args.target_key,
        # AnnData serializes mapping keys as HDF5 group names. Cell-type labels
        # can legitimately contain '/', so store the mapping as two aligned
        # arrays instead of using labels as dictionary keys.
        "target_labels": list(label_map),
        "target_codes": list(label_map.values()),
        "target_weights": args.weights, "neighbors": args.neighbors,
        "min_dist": args.min_dist, "seed": args.seed, "guided_cells": int((labels >= 0).sum()),
        "unlabeled_cells": int((labels < 0).sum()),
    }
    embedding_path = output / "methylvi_supervised_umap.h5ad"
    embedding.write_h5ad(embedding_path, compression="gzip")
    summary = {**embedding.uns["supervised_umap"], "target_mapping": label_map,
               "cells": embedding.n_obs, "embedding_h5ad": str(embedding_path),
               "coordinate_files": coordinate_files, "figure_files": figure_files}
    (output / "supervised_umap_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
