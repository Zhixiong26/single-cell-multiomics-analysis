#!/usr/bin/env python3
"""Adapt generated project configuration to the generic Scanpy implementation."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    sys.path.insert(0, str(root / "tools"))
    from _common import load_project, load_samples, write_json, write_tsv, SAMPLE_COLUMNS

    cfg = load_project(root)
    samples = [row for row in load_samples(root) if row["include"] and row["rna_path"]]
    if not samples:
        raise RuntimeError("no included RNA sample")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    annotation = cfg.get("annotation") or {}
    profile = annotation.get("profile")
    if profile and not Path(str(profile)).is_absolute():
        profile = str((root / str(profile)).resolve())
    resolved = {
        "project": {
            "organism": cfg["project"].get("organism", "unspecified"),
            "genome_build": cfg["project"].get("genome_id", "unspecified"),
        },
        "species": {"mitochondrial_prefixes": cfg["project"].get("mitochondrial_prefixes", ["MT-", "mt-"])},
        "analysis": {"scanpy": cfg["analysis"].get("scanpy", {})},
        "annotation": {"profile": profile, "unassigned_label": annotation.get("unassigned_label", "Unassigned")},
    }
    resolved_path = args.output_dir / "resolved_scanpy_config.json"
    samples_path = args.output_dir / "input_manifest.tsv"
    write_json(resolved_path, resolved)
    write_tsv(samples_path, SAMPLE_COLUMNS, samples)
    return subprocess.call([
        sys.executable, str(root / "Scripts" / "Scanpy" / "run_scanpy.py"),
        "--config", str(resolved_path), "--samples", str(samples_path),
        "--output-dir", str(args.output_dir),
    ], cwd=str(root))


if __name__ == "__main__":
    raise SystemExit(main())
