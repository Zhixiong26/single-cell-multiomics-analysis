#!/usr/bin/env python3
"""Export single-panel UMAPs from a completed VMR-MethylVI embedding."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(os.environ["VMR_MVI_RESULTS"]) / "methylvi_embedding.h5ad")
    parser.add_argument("--output-dir", type=Path, default=Path(os.environ["VMR_MVI_RESULTS"]))
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    embedding = ad.read_h5ad(args.input)
    if "X_umap" not in embedding.obsm:
        raise KeyError("VMR-MethylVI embedding lacks obsm['X_umap']")
    plots = (
        ("cell_type", "methylvi_vmr_umap_cell_type.png"),
        ("sample_id", "methylvi_vmr_umap_sample_id.png"),
        ("condition", "methylvi_vmr_umap_condition.png"),
        ("methylVI_leiden", "methylvi_vmr_umap_methylVI_leiden.png"),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for column, filename in plots:
        if column not in embedding.obs:
            raise KeyError(f"VMR-MethylVI embedding lacks obs[{column!r}]")
        figure = sc.pl.umap(
            embedding, color=column, show=False, return_fig=True,
            title=f"VMR-MethylVI UMAP — {column}",
        )
        figure.savefig(args.output_dir / filename, dpi=300, bbox_inches="tight")
    print(f"Wrote {len(plots)} VMR-MethylVI UMAPs to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
