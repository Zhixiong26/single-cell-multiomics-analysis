#!/usr/bin/env python3
"""Build a Scanpy VMR embedding from MethSCAn mean shrunken residuals."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.decomposition import PCA


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cell-metadata", type=Path)
    parser.add_argument("--expected-cell-ids", type=Path)
    parser.add_argument("--min-region-cell-fraction", type=float, default=0.05)
    parser.add_argument("--min-cell-regions", type=int, default=100)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.8)
    parser.add_argument("--impute-iterations", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def iterative_pca_impute(values, n_components, iterations, seed):
    missing = ~np.isfinite(values)
    col_means = np.nanmean(values, axis=0)
    col_means[~np.isfinite(col_means)] = 0.0
    filled = np.where(missing, col_means[None, :], values).astype(np.float32)
    mse = []
    if not missing.any():
        return filled, mse
    for _ in range(iterations):
        previous = filled[missing].copy()
        model = PCA(n_components=n_components, svd_solver="randomized", random_state=seed)
        scores = model.fit_transform(filled)
        reconstructed = model.inverse_transform(scores).astype(np.float32)
        filled[missing] = reconstructed[missing]
        mse.append(float(np.mean((filled[missing] - previous) ** 2)))
        if len(mse) > 1 and mse[-1] <= max(mse[0] * 0.001, 1e-10):
            break
    return filled, mse


def cell_metadata(cell_ids):
    samples, barcodes = [], []
    for cell_id in cell_ids:
        sample, separator, barcode = cell_id.partition("_")
        samples.append(sample if separator else "NA")
        barcodes.append(barcode if separator else cell_id)
    return pd.DataFrame(
        {"sample_id": samples, "barcode": barcodes}, index=pd.Index(cell_ids, name="cell_id")
    )


def main():
    args = parse_args()
    if not 0 < args.min_region_cell_fraction <= 1:
        raise ValueError("min-region-cell-fraction must be in (0, 1]")
    if args.min_cell_regions < 1 or args.n_pcs < 2 or args.neighbors < 2:
        raise ValueError("min-cell-regions must be positive; n-pcs/neighbors must be >= 2")
    if args.impute_iterations < 1 or args.leiden_resolution <= 0:
        raise ValueError("impute-iterations and leiden-resolution must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Scanpy output is not empty: %s" % args.output_dir)
    for child in ("tables", "figures", "objects"):
        (args.output_dir / child).mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.matrix, index_col=0)
    frame.index = frame.index.astype(str)

    # Cell IDs entering Scanpy must be exactly the canonical MethSCAn IDs
    # retained by methscan filter, in the same order.
    if args.expected_cell_ids:
        with args.expected_cell_ids.open() as handle:
            expected_ids = [line.strip() for line in handle if line.strip()]
        matrix_ids = frame.index.tolist()
        if matrix_ids != expected_ids:
            missing = sorted(set(expected_ids).difference(matrix_ids))
            extra = sorted(set(matrix_ids).difference(expected_ids))
            raise ValueError(
                "Matrix cell IDs do not exactly match filtered column_header.txt; "
                "missing=%s extra=%s"
                % (missing[:3], extra[:3])
            )
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError("Duplicate cell or VMR identifiers in matrix")
    values = frame.to_numpy(dtype=np.float32)
    input_finite = np.isfinite(values)
    region_observed_fraction = np.mean(input_finite, axis=0)
    region_keep = region_observed_fraction >= args.min_region_cell_fraction
    region_qc = pd.DataFrame({
        "observed_cell_fraction_before_cell_filter": region_observed_fraction,
        "pass_observed_cell_fraction": region_keep,
    }, index=pd.Index(frame.columns.astype(str), name="vmr"))
    values = values[:, region_keep]
    regions = frame.columns[region_keep].astype(str)
    cell_observed = np.sum(np.isfinite(values), axis=1)
    cell_keep = cell_observed >= args.min_cell_regions
    cell_qc = cell_metadata(frame.index.astype(str))
    cell_qc["n_observed_vmrs_after_region_filter"] = cell_observed
    cell_qc["observed_vmr_fraction_after_region_filter"] = (
        cell_observed / max(values.shape[1], 1)
    )
    cell_qc["pass_min_observed_vmrs"] = cell_keep
    values = values[cell_keep, :]
    cells = frame.index[cell_keep]
    if values.shape[0] < 3 or values.shape[1] < 2:
        raise ValueError("Too few cells or VMRs after missingness filters")
    n_pcs = min(args.n_pcs, values.shape[0] - 1, values.shape[1] - 1)
    if n_pcs < 2:
        raise ValueError("At least two PCA components are required")
    filled, imputation_mse = iterative_pca_impute(
        values, n_pcs, args.impute_iterations, args.seed
    )
    obs = cell_metadata(cells)
    obs["n_observed_vmrs"] = cell_observed[cell_keep]
    if args.cell_metadata:
        extra = pd.read_csv(args.cell_metadata, sep="\t", dtype=str, keep_default_na=False)
        if "cell_id" not in extra.columns or extra["cell_id"].duplicated().any():
            raise ValueError("cell metadata needs a unique cell_id column")
        extra = extra.set_index("cell_id")

        missing_metadata = obs.index.difference(extra.index)
        if len(missing_metadata):
            raise ValueError(
                "Cell metadata is missing retained cells, first: %s"
                % missing_metadata[0]
            )

        # Validate that the canonical cell_id agrees with the explicit sample
        # and barcode columns before joining annotation.
        if {"sample_id", "barcode"}.issubset(extra.columns):
            meta = extra.loc[obs.index, ["sample_id", "barcode"]]
            expected_ids = (
                meta["sample_id"].astype(str)
                + "_"
                + meta["barcode"].astype(str)
            )
            mismatch = expected_ids.index[expected_ids.to_numpy() != obs.index.to_numpy()]
            if len(mismatch):
                first = mismatch[0]
                raise ValueError(
                    "Metadata sample_id/barcode disagree with cell_id: %s"
                    % first
                )

        extra = extra.drop(
            columns=[c for c in ("sample_id", "barcode") if c in extra]
        )
        obs = obs.join(extra, how="left")
    var = pd.DataFrame(index=pd.Index(regions, name="vmr"))
    var["observed_cell_fraction_after_cell_filter"] = np.mean(np.isfinite(values), axis=0)
    adata = AnnData(X=filled, obs=obs, var=var)
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=args.seed)
    adata.obsm["X_pca"] = pca.fit_transform(filled).astype(np.float32)
    adata.uns["pca_variance_ratio"] = pca.explained_variance_ratio_.tolist()
    sc.pp.neighbors(adata, n_neighbors=min(args.neighbors, adata.n_obs - 1), use_rep="X_pca", random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.leiden_resolution, key_added="leiden", random_state=args.seed)

    coordinates = adata.obs.copy()
    coordinates["PC1"] = adata.obsm["X_pca"][:, 0]
    coordinates["PC2"] = adata.obsm["X_pca"][:, 1]
    coordinates["UMAP1"] = adata.obsm["X_umap"][:, 0]
    coordinates["UMAP2"] = adata.obsm["X_umap"][:, 1]
    coordinates.to_csv(args.output_dir / "tables" / "cell_embedding.tsv", sep="\t")
    cell_qc.to_csv(args.output_dir / "tables" / "cell_missingness_qc.tsv", sep="\t")
    region_qc.to_csv(args.output_dir / "tables" / "vmr_missingness_qc.tsv", sep="\t")
    adata.var.to_csv(args.output_dir / "tables" / "vmr_qc.tsv", sep="\t")
    summary = (
        adata.obs.groupby(["sample_id", "leiden"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    summary.to_csv(args.output_dir / "tables" / "sample_cluster_counts.tsv", sep="\t", index=False)

    plot_colors = ["sample_id", "leiden"]
    if "rna_cell_type" in adata.obs:
        plot_colors.append("rna_cell_type")
    for color in plot_colors:
        fig, ax = plt.subplots(figsize=(7, 6))
        category_values = adata.obs[color].astype(str)
        # Matplotlib's default color cycle has only 10 colors, which caused
        # distinct cell types to be rendered with identical colors (the
        # current annotation contains 14 categories).  Use a deterministic
        # discrete palette with enough visually distinct colors instead.
        categories = sorted(category_values.unique())
        palette = list(plt.get_cmap("tab20").colors)
        for i, category in enumerate(categories):
            mask = category_values == category
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=5,
                alpha=0.7,
                color=palette[i % len(palette)],
                label=category,
            )
        ax.set(xlabel="UMAP1", ylabel="UMAP2", title="VMR methylation: %s" % color)
        ax.legend(markerscale=2, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        fig.tight_layout()
        fig.savefig(args.output_dir / "figures" / ("umap_%s.png" % color), dpi=300)
        plt.close(fig)
    adata.write_h5ad(args.output_dir / "objects" / "methscan_vmr.h5ad", compression="gzip")
    parameters = vars(args).copy()
    parameters.update({
        "input_cells": int(frame.shape[0]), "input_regions": int(frame.shape[1]),
        "retained_cells": int(adata.n_obs), "retained_regions": int(adata.n_vars),
        "effective_n_pcs": int(n_pcs), "imputation_mse": imputation_mse,
        "input_missing_fraction": float(1 - np.mean(input_finite)),
        "retained_matrix_missing_fraction_before_imputation": float(
            1 - np.mean(np.isfinite(values))
        ),
        # Labels come from the pre-filter Scanpy whitelist; no second cell-type
        # selection is performed after MethSCAn filtering.
        "rna_annotation_available_cells": int(
            adata.obs["rna_cell_type"].notna().sum()
        ) if "rna_cell_type" in adata.obs else None,
        "scanpy_version": sc.__version__,
    })
    for key, value in list(parameters.items()):
        if isinstance(value, Path):
            parameters[key] = str(value.resolve())
    with (args.output_dir / "run_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
