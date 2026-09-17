#!/usr/bin/env python3
"""Create a nested top-N MethylVI H5MU from a ranked integer-count input."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mudata
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--features", type=int, required=True)
    args = parser.parse_args()
    if args.features < 2 or not args.input.is_file():
        raise ValueError("A valid input and at least two features are required")
    source = mudata.read_h5mu(args.input)
    adata = source["mCG"]
    if not {"mc", "cov"}.issubset(adata.layers):
        raise KeyError("mCG input requires integer mc/cov layers")
    if args.features > adata.n_vars:
        raise ValueError(f"Requested {args.features:,} from only {adata.n_vars:,} features")
    if "selection_rank" not in adata.var:
        raise KeyError("mCG var requires deterministic selection_rank")
    ranks = np.asarray(adata.var["selection_rank"], dtype=np.int64)
    order = np.argsort(ranks, kind="stable")
    selected = order[: args.features]
    if set(ranks[selected]) != set(range(1, args.features + 1)):
        raise ValueError("selection_rank is not a complete nested top-N ranking")
    subset = adata[:, selected].copy()
    if subset.n_vars != args.features:
        raise RuntimeError("Nested feature hard check failed")
    output = mudata.MuData({"mCG": subset})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp.h5mu")
    temporary.unlink(missing_ok=True)
    output.write_h5mu(temporary, compression="gzip")
    temporary.replace(args.output)
    summary = {
        "source": str(args.input.resolve()),
        "output": str(args.output.resolve()),
        "cells": int(subset.n_obs),
        "source_features": int(adata.n_vars),
        "selected_features": int(subset.n_vars),
        "selection": "selection_rank <= requested features",
    }
    args.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
