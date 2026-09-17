#!/usr/bin/env python3
"""Run generic multi-sample 10x RNA QC, integration, clustering, and guarded annotation."""

from __future__ import print_function

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Scripts.Common.annotation import evaluate_annotation_profile


def load_table(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_structured(path):
    text = Path(path).read_text(encoding="utf-8")
    try:
        import yaml
    except ImportError:
        return json.loads(text)
    return yaml.safe_load(text)


def extract_matrix(path, root, sample_id, rna_format):
    path = Path(path)
    if rna_format == "10x_h5":
        if not path.is_file():
            raise ValueError("10x H5 input is not a file: %s" % path)
        return path
    if path.is_dir():
        return path
    if not zipfile.is_zipfile(str(path)):
        raise ValueError("RNA input must be a 10x directory or ZIP: %s" % path)
    target = root / sample_id
    if target.exists():
        shutil.rmtree(str(target))
    target.mkdir(parents=True)
    with zipfile.ZipFile(str(path)) as archive:
        archive.extractall(str(target))
    candidates = []
    for marker in ("matrix.mtx", "matrix.mtx.gz"):
        candidates.extend(item.parent for item in target.rglob(marker))
    unique = sorted(set(candidates))
    if len(unique) != 1:
        raise ValueError("expected exactly one 10x matrix in %s, found %d" % (path, len(unique)))
    return unique[0]


def analysis_signature(cfg, samples):
    payload = {
        "organism": cfg["project"]["organism"],
        "genome_build": cfg["project"]["genome_build"],
        "scanpy": cfg.get("analysis", {}).get("scanpy", {}),
        "samples": [{key: row.get(key, "") for key in ("sample_id", "batch", "rna_path", "cell_id_prefix")} for row in samples if row.get("rna_path")],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def apply_annotation(adata, cfg, signature, output_dir):
    label = cfg.get("annotation", {}).get("unassigned_label", "Unassigned")
    observed = sorted(adata.obs["leiden"].astype(str).unique())
    profile_path = cfg.get("annotation", {}).get("profile")
    profile = load_structured(profile_path) if profile_path else None
    annotations, status, template = evaluate_annotation_profile(observed, signature, profile, label)
    if annotations:
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(lambda item: annotations[item]["cell_type"])
        adata.obs["annotation_confidence"] = adata.obs["leiden"].astype(str).map(lambda item: annotations[item].get("confidence", "reviewed"))
    else:
        adata.obs["cell_type"] = label
        adata.obs["annotation_confidence"] = "unreviewed"
    (output_dir / "annotation_template.yaml").write_text(json.dumps(template, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "annotation_guard_status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return status


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--samples", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    samples = load_table(args.samples)
    rna_samples = [row for row in samples if row.get("rna_path")]
    if not rna_samples:
        raise SystemExit("input manifest contains no RNA samples")
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit("Scanpy stage environment is incomplete: %s" % exc)

    settings = cfg.get("analysis", {}).get("scanpy", {})
    species = cfg.get("species", {})
    qc = settings.get("qc", {})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    extracted = args.output_dir / "extracted_inputs"
    adatas = []
    qc_rows = []
    for row in rna_samples:
        matrix_path = extract_matrix(row["rna_path"], extracted, row["sample_id"], row.get("rna_format", "10x_mtx"))
        if row.get("rna_format") == "10x_h5":
            sample = sc.read_10x_h5(str(matrix_path))
            sample.var_names_make_unique()
        else:
            sample = sc.read_10x_mtx(str(matrix_path), var_names="gene_symbols", make_unique=True)
        prefix = row.get("cell_id_prefix") or row["sample_id"]
        sample.obs_names = [name if name.startswith(prefix + "_") else prefix + "_" + name for name in sample.obs_names]
        # The stable samples.tsv contract intentionally keeps only metadata
        # required by every project.  Do not assume donor/tissue columns exist.
        for key in ("sample_id", "condition", "batch"):
            sample.obs[key] = row.get(key, "")
        sample.layers["counts"] = sample.X.copy()
        upper = sample.var_names.str.upper()
        mt_prefixes = tuple(item.upper() for item in species.get("mitochondrial_prefixes", ["MT-"]))
        ribo_prefixes = tuple(item.upper() for item in species.get("ribosomal_prefixes", ["RPL", "RPS"]))
        sample.var["mt"] = upper.str.startswith(mt_prefixes)
        sample.var["ribo"] = upper.str.startswith(ribo_prefixes)
        sc.pp.calculate_qc_metrics(sample, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
        keep = (
            (sample.obs["n_genes_by_counts"] >= int(qc.get("min_genes", 200))) &
            (sample.obs["n_genes_by_counts"] < int(qc.get("max_genes", 6000))) &
            (sample.obs["total_counts"] >= int(qc.get("min_counts", 500))) &
            (sample.obs["pct_counts_mt"] < float(qc.get("max_mt_percent", 10.0)))
        )
        before = sample.n_obs
        sample = sample[keep].copy()
        predicted = np.zeros(sample.n_obs, dtype=bool)
        if settings.get("scrublet", {}).get("enabled", True) and sample.n_obs >= int(settings.get("scrublet", {}).get("min_cells", 100)):
            try:
                import scrublet as scr
                expected = float(settings.get("scrublet", {}).get("expected_doublet_rate", 0.06))
                scores, predicted = scr.Scrublet(sample.layers["counts"], expected_doublet_rate=expected).scrub_doublets()
                sample.obs["doublet_score"] = scores
            except Exception as exc:
                raise RuntimeError("Scrublet failed for %s: %s" % (row["sample_id"], exc))
        sample.obs["predicted_doublet"] = predicted
        sample = sample[~sample.obs["predicted_doublet"]].copy()
        sc.pp.filter_genes(sample, min_cells=int(qc.get("min_cells_per_gene", 3)))
        qc_rows.append({"sample_id": row["sample_id"], "input_cells": before, "retained_cells": sample.n_obs, "predicted_doublets": int(predicted.sum())})
        adatas.append(sample)

    combined = ad.concat(adatas, join="outer", merge="same", fill_value=0)
    sc.pp.normalize_total(combined, target_sum=float(settings.get("target_sum", 10000)))
    sc.pp.log1p(combined)
    combined.raw = combined
    n_hvg = min(int(settings.get("n_hvg", 2000)), combined.n_vars)
    sc.pp.highly_variable_genes(combined, n_top_genes=n_hvg, flavor=settings.get("hvg_flavor", "seurat"), batch_key=settings.get("batch_key", "batch"))
    work = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(int(settings.get("pca_components", 50)), work.n_obs - 1, work.n_vars - 1), random_state=int(settings.get("seed", 0)))
    representation = "X_pca"
    if settings.get("integration", "harmony") == "harmony":
        try:
            import scanpy.external as sce
            sce.pp.harmony_integrate(work, key=settings.get("batch_key", "batch"), basis="X_pca", adjusted_basis="X_pca_harmony")
        except Exception as exc:
            raise RuntimeError("Harmony integration failed: %s" % exc)
        representation = "X_pca_harmony"
    sc.pp.neighbors(work, n_neighbors=int(settings.get("n_neighbors", 15)), n_pcs=min(int(settings.get("n_pcs", 30)), work.obsm[representation].shape[1]), use_rep=representation, random_state=int(settings.get("seed", 0)))
    sc.tl.umap(work, random_state=int(settings.get("seed", 0)))
    sc.tl.leiden(work, resolution=float(settings.get("leiden_resolution", 0.8)), random_state=int(settings.get("seed", 0)))
    sc.tl.rank_genes_groups(work, "leiden", method="wilcoxon", use_raw=True)
    markers = sc.get.rank_genes_groups_df(work, group=None)
    markers.to_csv(args.output_dir / "ranked_markers.tsv.gz", sep="\t", index=False)
    signature = analysis_signature(cfg, rna_samples)
    annotation_status = apply_annotation(work, cfg, signature, args.output_dir)
    work.obs[["sample_id", "condition", "batch", "leiden", "cell_type", "annotation_confidence"]].to_csv(args.output_dir / "cell_id_cell_type.tsv", sep="\t", index_label="cell_id")
    pd.DataFrame(qc_rows).to_csv(args.output_dir / "sample_qc.tsv", sep="\t", index=False)
    work.write_h5ad(args.output_dir / "scanpy_result.h5ad", compression="gzip")
    summary = {
        "status": "complete", "stage": "scanpy", "analysis_signature": signature,
        "annotation_status": annotation_status["status"], "n_cells": work.n_obs,
        "n_genes": work.n_vars, "n_clusters": int(work.obs["leiden"].nunique()),
        "versions": {"python": sys.version.split()[0], "scanpy": sc.__version__, "anndata": ad.__version__},
    }
    (args.output_dir / "completion_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
