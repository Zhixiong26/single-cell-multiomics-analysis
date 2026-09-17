#!/usr/bin/env python3
"""Export named, single-panel pre-MethylVI ALLCools embeddings."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(os.environ["SCMO_ALLCOOLS_H5AD"]))
    parser.add_argument("--output-dir", type=Path, default=Path(os.environ["SCMO_ALLCOOLS_ROOT"]))
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    adata = ad.read_h5ad(args.input)
    for key in ("X_umap", "X_tsne"):
        if key not in adata.obsm:
            raise KeyError(f"ALLCools H5AD lacks obsm[{key!r}]")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for basis, column, filename in (
        ("tsne", "L1", "allcools_5kb_tsne_L1.png"),
        ("umap", "L1", "allcools_5kb_umap_L1.png"),
        ("umap", "manual_celltype", "allcools_original_embedding_cell_type.png"),
        ("umap", "cohort", "allcools_original_embedding_cohort.png"),
    ):
        if column not in adata.obs:
            raise KeyError(f"ALLCools H5AD lacks obs[{column!r}]")
        figure = sc.pl.embedding(
            adata, basis=basis, color=column, show=False, return_fig=True,
            title=f"ALLCools 5-kb {basis.upper()} — {column}",
        )
        figure.savefig(args.output_dir / filename, dpi=300, bbox_inches="tight")
    print(f"Wrote 4 single-panel ALLCools embeddings to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
