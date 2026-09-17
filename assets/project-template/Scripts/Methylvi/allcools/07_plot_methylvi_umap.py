#!/usr/bin/env python3
"""Export single-panel, post-MethylVI UMAPs from the finished embedding."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(os.environ["SCMO_MVI_RESULTS"]) / "methylvi_embedding.h5ad")
    parser.add_argument("--output-dir", type=Path, default=Path(os.environ["SCMO_MVI_RESULTS"]))
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    adata = ad.read_h5ad(args.input)
    if "X_umap" not in adata.obsm:
        raise KeyError("MethylVI embedding lacks obsm['X_umap']")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plots = (
        ("cell_type", "cell_type"),
        ("sample_id", "sample_id"),
        ("condition", "condition"),
        ("methylVI_leiden", "methylVI_leiden"),
    )
    for column, output_name in plots:
        if column not in adata.obs:
            raise KeyError(f"MethylVI embedding lacks obs[{column!r}]")
        figure = sc.pl.umap(
            adata, color=column, show=False, return_fig=True,
            title=f"MethylVI UMAP — {column}",
        )
        figure.savefig(args.output_dir / f"methylvi_umap_{output_name}.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {len(plots)} single-panel UMAPs to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
