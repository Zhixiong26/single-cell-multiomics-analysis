#!/usr/bin/env python3
"""Validate and summarize all VMR plus pooled-DMR MethylVI models."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--expected-cells", type=int, required=True)
    parser.add_argument("--threshold", action="append", required=True)
    parser.add_argument("--feature-target", action="append", type=int, required=True)
    args = parser.parse_args()
    dmr_summary = json.loads((args.root / "shared/dmr/prepare_summary.json").read_text())
    if dmr_summary.get("status") != "complete" or dmr_summary.get("top_n_truncation") is not None:
        raise ValueError("Pooled DMR preparation is incomplete or unexpectedly Top-N truncated")
    models = []
    for threshold in args.threshold:
        for base_features in args.feature_target:
            route = args.root / "vmr_dmr" / f"var_{threshold}" / f"features_{base_features}"
            input_summary = json.loads((route / "input.summary.json").read_text())
            run_summary = json.loads((route / "results/run_summary.json").read_text())
            required = [
                route / "input.COMPLETE", route / "model.COMPLETE",
                route / "results/methylvi_embedding.h5ad",
                route / "results/methylvi_vmr_umap_cell_type.png",
                route / "results/methylvi_vmr_umap_sample_id.png",
                route / "results/methylvi_vmr_umap_condition.png",
                route / "results/methylvi_vmr_umap_methylVI_leiden.png",
            ]
            missing = [str(path) for path in required if not path.is_file()]
            if missing:
                raise FileNotFoundError("Missing model artifacts: " + ", ".join(missing))
            embedding = ad.read_h5ad(route / "results/methylvi_embedding.h5ad", backed="r")
            if embedding.n_obs != args.expected_cells or embedding.n_vars != 20:
                raise ValueError(f"Unexpected embedding shape for {route}: {embedding.shape}")
            if input_summary["combined_features"] != run_summary["features"]:
                raise ValueError(f"Input/training feature mismatch for {route}")
            models.append({
                "variance_threshold": float(threshold), "base_vmr_features": base_features,
                **{key: input_summary[key] for key in (
                    "nonoverlapping_dmrs_added", "dmrs_overlapping_vmr_removed", "combined_features"
                )},
                "embedding": str((route / "results/methylvi_embedding.h5ad").resolve()),
                "status": "complete",
            })
    summary = {
        "status": "complete", "cells": args.expected_cells,
        "dmr_selection": "all exact-unique hypo-DMRs with raw p < 0.01 and abs methdiff >= 0.25; overlap-merged; no Top-N truncation",
        "vmr_dmr_overlap_rule": "retain selected VMRs and append only zero-overlap pooled DMRs",
        "models": models,
    }
    (args.root / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
