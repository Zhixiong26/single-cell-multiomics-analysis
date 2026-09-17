#!/usr/bin/env python3
"""Train MethylVI on raw counts for ALLCools-retained 5-kb bins."""
from __future__ import annotations

import argparse
import json
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import mudata
import numpy as np
import scanpy as sc
import scvi
import torch
from scvi.external import METHYLVI


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(os.environ["SCMO_MVI_INPUT"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["SCMO_MVI_RESULTS"]))
    parser.add_argument("--batch-key", default=os.environ["SCMO_BATCH_KEY"])
    parser.add_argument("--epochs", type=int, default=int(os.environ["SCMO_EPOCHS"]))
    parser.add_argument("--batch-size", type=int, default=int(os.environ["SCMO_BATCH_SIZE"]))
    parser.add_argument("--threads", type=int, default=int(os.environ["SCMO_THREADS"]))
    parser.add_argument("--seed", type=int, default=int(os.environ["SCMO_SEED"]))
    parser.add_argument("--accelerator", choices=("auto", "cpu", "gpu"), default="auto")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = args.seed
    torch.set_num_threads(args.threads)

    mdata = mudata.read_h5mu(args.input)
    adata = mdata["mCG"]
    if args.batch_key not in adata.obs:
        raise KeyError(f"Batch key absent: {args.batch_key}")
    if adata.obs[args.batch_key].nunique() < 2:
        warnings.warn(f"Batch key {args.batch_key!r} has fewer than two levels")
    mc, cov = adata.layers["mc"], adata.layers["cov"]
    if not np.issubdtype(mc.dtype, np.integer) or not np.issubdtype(cov.dtype, np.integer):
        raise ValueError("mc/cov layers must contain integer counts")
    if np.any(mc > cov):
        raise ValueError("Invalid count layers: mc > cov")

    METHYLVI.setup_mudata(
        mdata, mc_layer="mc", cov_layer="cov", batch_key=args.batch_key,
        methylation_contexts=["mCG"], modalities={"batch_key": "mCG"},
    )
    model = METHYLVI(
        mdata, n_latent=20, n_hidden=128, n_layers=1,
        likelihood="betabinomial", dispersion="region",
    )
    accelerator = args.accelerator
    if accelerator == "auto":
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    model.train(
        max_epochs=args.epochs, early_stopping=True, batch_size=args.batch_size,
        accelerator=accelerator, devices=1,
    )
    model.save(args.output / "model", overwrite=True, save_anndata=False)
    latent = np.asarray(model.get_latent_representation(batch_size=args.batch_size), dtype=np.float32)
    np.save(args.output / "methylvi_latent.npy", latent)
    embedding = ad.AnnData(X=latent, obs=adata.obs.copy())
    embedding.obsm["X_methylVI"] = latent
    sc.pp.neighbors(embedding, use_rep="X_methylVI", n_neighbors=15, random_state=args.seed)
    sc.tl.umap(embedding, random_state=args.seed)
    sc.tl.leiden(
        embedding, resolution=1.0, key_added="methylVI_leiden", random_state=args.seed
    )
    embedding.write_h5ad(args.output / "methylvi_embedding.h5ad", compression="gzip")
    summary = {
        "input": str(args.input), "cells": embedding.n_obs, "features": adata.n_vars,
        "batch_key": args.batch_key, "batch_levels": int(adata.obs[args.batch_key].nunique()),
        "accelerator": accelerator, "epochs_requested": args.epochs,
    }
    (args.output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
