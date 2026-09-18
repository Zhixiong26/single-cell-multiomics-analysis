#!/usr/bin/env python3
"""Execute a built-in workflow task from generated configuration."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--task", required=True)
    parser.add_argument("--task-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    sys.path.insert(0, str(root / "tools"))
    from _common import WorkflowError, load_environments, load_project, project_files, resolve_path, write_json

    cfg = load_project(root)
    plan = json.loads((root / ".workflow" / "runs" / args.run_id / "plan.json").read_text())
    item = next(value for value in plan["tasks"] if value["id"] == args.task)
    result = Path(plan["result_dir"])
    result.mkdir(parents=True, exist_ok=True)
    args.task_dir.mkdir(parents=True, exist_ok=True)
    files = project_files(root)
    analysis = cfg["analysis"]
    meth = analysis.get("methscan", {})
    dmr = meth.get("dmr", {})
    allcools = analysis.get("allcools", {})
    mvi = analysis.get("methylvi", {})
    annotation = cfg.get("annotation") or {}
    annotation_table = resolve_path(root, annotation.get("table") or annotation.get("path"))
    generated_annotation = result / "annotation.tsv"
    if annotation_table is None and generated_annotation.is_file():
        annotation_table = generated_annotation
    if annotation_table and annotation_table.is_file() and (
            annotation.get("cell_id_column", "cell_id") != "cell_id"
            or annotation.get("cell_type_column", "cell_type") != "cell_type"):
        normalized = result / "annotation.normalized.tsv"
        with annotation_table.open(newline="") as source, normalized.open("w", newline="") as target:
            rows = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(target, fieldnames=["cell_id", "cell_type"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({"cell_id": row[annotation.get("cell_id_column", "cell_id")],
                                 "cell_type": row[annotation.get("cell_type_column", "cell_type")]})
        annotation_table = normalized
    references = cfg.get("references") or {}
    chrom_sizes = resolve_path(root, references.get("chrom_sizes"))
    blacklist = resolve_path(root, references.get("blacklist"))
    if args.task != "scanpy" and (not chrom_sizes or not blacklist):
        raise WorkflowError("built-in methylation tasks require chrom_sizes and blacklist")

    environments = {row["stage"]: row for row in load_environments(root)}
    def python_for(stage):
        value = environments.get(stage, {}).get("python", "")
        if not value:
            raise WorkflowError("environment Python is absent for %s" % stage)
        return str(resolve_path(root, value))

    def executable_for(stage):
        value = environments.get(stage, {}).get("executable", "")
        if not value:
            raise WorkflowError("environment executable is absent for %s" % stage)
        return str(resolve_path(root, value))

    scanpy_python = python_for("scanpy_allcools")
    methscan_python = python_for("methscan") if args.task != "scanpy" else scanpy_python
    methylvi_python = python_for("methylvi") if args.task != "scanpy" else scanpy_python
    methscan_exe = executable_for("methscan") if args.task != "scanpy" else scanpy_python
    allcools_prefix = str(Path(scanpy_python).resolve().parent.parent)
    methylvi_prefix = str(Path(methylvi_python).resolve().parent.parent)
    meth_root = result / "methscan"
    mvi_root = result / "methylvi"
    feature_targets = [int(value) for value in mvi.get("feature_targets", [10000, 30000])]
    max_features = max(feature_targets)
    thresholds = [str(value) for value in meth.get("vmr_thresholds", [0.01, 0.02, 0.05])]

    env = os.environ.copy()
    env.update({
        "SCMO_PROJECT_ROOT": str(root), "SCMO_RESULT_DIR": str(result),
        "SCMO_MANAGED_RUN_ID": args.run_id,
        "SCMO_METHSCAN_RUN_DIR": str(meth_root), "SCMO_SAMPLES_TSV": str(files["samples"]),
        "SCMO_ANNOTATION": str(annotation_table) if annotation_table else "", "SCMO_CHROM_SIZES": str(chrom_sizes or ""),
        "SCMO_BLACKLIST": str(blacklist), "SCMO_METHSCAN_EXE": methscan_exe,
        "SCMO_METHSCAN_PYTHON": methscan_python, "SCMO_SCANPY_PYTHON": scanpy_python,
        "SCMO_ALLCOOLS_ENV": allcools_prefix, "SCMO_METHYLVI_ENV": methylvi_prefix,
        "SCMO_ALLC_SOURCE": next((row["allc_root"] for row in __import__("_common").load_samples(root) if row["include"] and row["allc_root"]), str(root)),
        "SCMO_SAMPLE_IDS": " ".join(row["sample_id"] for row in __import__("_common").load_samples(root) if row["include"] and row["allc_root"]),
        "SCMO_ANNOTATION_APPROVED": "1" if annotation.get("review_status") == "approved" else "0",
        "SCMO_MIN_SITES": str(meth.get("min_sites", 300000)),
        "SCMO_MIN_METH_PERCENT": str(meth.get("min_meth_percent", 50)),
        "SCMO_MAX_METH_PERCENT": str(meth.get("max_meth_percent", 100)),
        "SCMO_SMOOTH_BANDWIDTH": str(meth.get("smooth_bandwidth", 1000)),
        "SCMO_SCAN_BANDWIDTH": str(meth.get("scan_bandwidth", 2000)),
        "SCMO_SCAN_STEPSIZE": str(meth.get("scan_stepsize", 100)),
        "SCMO_VMR_THRESHOLDS": " ".join(thresholds), "SCMO_MIN_CELLS": str(meth.get("min_cells", 6)),
        "SCMO_HYPO_RAW_P": str(dmr.get("raw_p", 0.01)),
        "SCMO_HYPO_MIN_ABS_DIFF": str(dmr.get("min_abs_diff", 0.25)),
        "SCMO_HYPO_TOP_PER_CELL_TYPE": str(dmr.get("top_per_cell_type", 200)),
        "SCMO_MC_CONTEXT": str(allcools.get("mc_context", "CGN")),
        "SCMO_BIN_SIZE": str(allcools.get("bin_size", 5000)),
        "SCMO_BLACKLIST_FRACTION": str(allcools.get("blacklist_fraction", 0.2)),
        "SCMO_FEATURE_TARGETS": " ".join(map(str, feature_targets)),
        "SCMO_TARGET_FEATURES": str(max_features), "SCMO_MAX_FEATURES": str(max_features),
        "SCMO_EPOCHS": str(mvi.get("epochs", 500)),
        "SCMO_BATCH_SIZE": str(mvi.get("batch_size", 32)), "SCMO_SEED": str(mvi.get("seed", 0)),
        "SCMO_THREADS": os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("SCMO_CPUS", "1")),
        "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    })

    def run(command, extra=None):
        task_env = dict(env)
        task_env.update(extra or {})
        subprocess.check_call([str(value) for value in command], cwd=str(root), env=task_env)

    def count_cells():
        path = meth_root / "03_filtered" / "column_header.txt"
        if not path.is_file():
            raise WorkflowError("filtered MethSCAn cell list is absent: %s" % path)
        count = sum(bool(line.strip()) for line in path.open())
        if count < 1:
            raise WorkflowError("filtered MethSCAn cell list is empty")
        return count

    def record(paths):
        resolved = [str(Path(path).resolve()) for path in paths]
        missing = [path for path in resolved if not Path(path).exists()]
        if missing:
            raise WorkflowError("task outputs are absent: %s" % ", ".join(missing))
        write_json(args.task_dir / "task_outputs.json", {"task": args.task, "artifacts": resolved})

    scripts = root / "Scripts"
    task = args.task
    parameters = item.get("parameters", {})
    if task == "scanpy":
        out = result / "scanpy"
        run([scanpy_python, scripts / "Scanpy/run_scanpy_project.py", "--project", root, "--output-dir", out])
        record([out / "completion_summary.json", out / "scanpy_result.h5ad"])
    elif task == "methscan_select_convert":
        run(["bash", scripts / "Methscan/run_methscan_common.sbatch", meth_root])
        if annotation_table is None:
            with (meth_root / "00_scanpy_selected/input_manifest.tsv").open(newline="") as handle:
                selected = list(csv.DictReader(handle, delimiter="\t"))
            with generated_annotation.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["cell_id", "cell_type"], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in selected:
                    writer.writerow({"cell_id": row["cell_id"], "cell_type": row.get("rna_cell_type") or "Unassigned"})
        record([meth_root / "common.COMPLETE", meth_root / "01_cov/input_manifest.tsv"])
    elif task in {"methscan_prepare", "methscan_filter", "methscan_smooth"}:
        stage = task.split("_", 1)[1]
        run(["bash", scripts / "Methscan/run_methscan_qc_stage.sbatch", meth_root, stage])
        record([meth_root / (stage + ".COMPLETE")])
    elif task.startswith("methscan_vmr_"):
        threshold = str(parameters["threshold"])
        run(["bash", scripts / "Methscan/run_methscan_branch.sbatch", meth_root, threshold])
        record([meth_root / ("var_%s.COMPLETE" % threshold), meth_root / "04_scan" / ("var_%s" % threshold) / "VMRs.bed"])
    elif task == "methscan_pairwise_dmr":
        run(["bash", scripts / "Methscan/run_methscan_methdiff.sbatch", meth_root])
        record([meth_root / "methdiff.COMPLETE", meth_root / "07_methdiff/pairwise_summary.tsv"])
    elif task == "methscan_hypo_heatmaps":
        run(["bash", scripts / "Methscan/run_methscan_hypo_heatmaps.sbatch", meth_root])
        record([meth_root / "hypo_heatmaps.COMPLETE"])
    elif task == "methscan_pooled_dmr":
        run(["bash", scripts / "Methscan/run_methscan_pooled_methdiff.sbatch", meth_root])
        record([meth_root / "pooled/methdiff.COMPLETE", meth_root / "pooled/07_methdiff/pairwise_summary.tsv"])
    elif task == "allcools_features":
        route, shared = mvi_root / "allcools", mvi_root / "allcools/shared"
        maximum = route / ("features_%d" % max_features)
        expected = count_cells()
        extra = {
            "SCMO_METHSCAN_RUN_DIR": str(meth_root), "SCMO_INPUT_MANIFEST": str(meth_root / "00_scanpy_selected/input_manifest.tsv"),
            "SCMO_FILTERED_CELL_IDS": str(meth_root / "03_filtered/column_header.txt"),
            "SCMO_COV_DIR": str(meth_root / "01_cov/cov"), "SCMO_EXISTING_ALLC_DIR": str(shared / "input_allc"),
            "SCMO_EXPECTED_CELLS": str(expected), "SCMO_MVI_ROOT": str(shared),
            "SCMO_ALLCOOLS_ROOT": str(shared / "allcools_features"), "SCMO_ALLC_DIR": str(shared / "input_allc"),
            "SCMO_ALLC_TABLE": str(shared / "selected_cells.allc.tsv"), "SCMO_MCDS": str(shared / "mc.mcds"),
            "SCMO_ALLCOOLS_H5AD": str(maximum / "allcools.h5ad"), "SCMO_MVI_INPUT": str(maximum / "input.h5mu"),
            "SCMO_MVI_RESULTS": str(maximum / "results"),
        }
        run(["bash", scripts / "Methylvi/allcools/run.sh", "prepare"], extra)
        run(["bash", scripts / "Methylvi/allcools/run.sh", "cluster"], extra)
        run(["bash", scripts / "Methylvi/allcools/run.sh", "build"], extra)
        for target in feature_targets:
            if target != max_features:
                run([methylvi_python, scripts / "Methylvi/shared/06_subset_methylvi_input.py", "--input", maximum / "input.h5mu",
                     "--output", route / ("features_%d/input.h5mu" % target), "--features", target], extra)
        record([route / ("features_%d/input.h5mu" % value) for value in feature_targets])
    elif task.startswith("methylvi_allcools_"):
        features = int(parameters["features"])
        route = mvi_root / "allcools" / ("features_%d" % features)
        train_model(root, scripts, methylvi_python, route, env, run)
        record([route / "results/methylvi_embedding.h5ad", route / "model.COMPLETE"])
    elif task.startswith("methylvi_vmr_features_"):
        threshold = str(parameters["threshold"])
        route, shared = mvi_root / "vmr" / ("var_%s" % threshold), mvi_root / "vmr" / ("var_%s/shared" % threshold)
        maximum = route / ("features_%d" % max_features)
        extra = vmr_environment(env, root, meth_root, annotation_table, chrom_sizes, blacklist,
                                methylvi_prefix, threshold, shared, maximum, max_features, count_cells(), allcools.get("mc_context", "CGN"))
        run(["bash", scripts / "Methylvi/vmr/run.sh", "prepare"], extra)
        run(["bash", scripts / "Methylvi/vmr/run.sh", "build"], extra)
        for target in feature_targets:
            if target != max_features:
                run([methylvi_python, scripts / "Methylvi/shared/06_subset_methylvi_input.py", "--input", maximum / "input.h5mu",
                     "--output", route / ("features_%d/input.h5mu" % target), "--features", target], extra)
        record([route / ("features_%d/input.h5mu" % value) for value in feature_targets])
    elif task.startswith("methylvi_vmr_dmr_"):
        threshold, features = str(parameters["threshold"]), int(parameters["features"])
        route = mvi_root / "vmr_dmr" / ("var_%s" % threshold) / ("features_%d" % features)
        base = mvi_root / "vmr" / ("var_%s" % threshold) / ("features_%d/input.h5mu" % features)
        dmr_input = mvi_root / "vmr_dmr/shared/dmr/input.h5mu"
        run([methylvi_python, scripts / "Methylvi/vmr_dmr/02_combine_vmr_and_dmr_inputs.py",
             "--vmr-input", base, "--dmr-input", dmr_input, "--output", route / "input.h5mu"])
        train_model(root, scripts, methylvi_python, route, env, run)
        record([route / "input.h5mu", route / "results/methylvi_embedding.h5ad", route / "model.COMPLETE"])
    elif task.startswith("methylvi_vmr_"):
        threshold, features = str(parameters["threshold"]), int(parameters["features"])
        route = mvi_root / "vmr" / ("var_%s" % threshold) / ("features_%d" % features)
        train_model(root, scripts, methylvi_python, route, env, run)
        record([route / "results/methylvi_embedding.h5ad", route / "model.COMPLETE"])
    elif task == "pooled_dmr_prepare":
        out = mvi_root / "vmr_dmr/shared/dmr"
        run([methylvi_python, scripts / "Methylvi/vmr_dmr/01_prepare_all_unique_pooled_dmrs.py",
             "--pairwise-summary", meth_root / "pooled/07_methdiff/pairwise_summary.tsv",
             "--blacklist", blacklist, "--chrom-sizes", chrom_sizes,
             "--blacklist-fraction", allcools.get("blacklist_fraction", 0.2),
             "--raw-p", dmr.get("raw_p", 0.01), "--min-abs-diff", dmr.get("min_abs_diff", 0.25),
             "--sort-threads", env["SCMO_THREADS"], "--output-dir", out])
        record([out / "all_unique_hypo_DMRs.merged.bed", out / "prepare_summary.json"])
    elif task == "pooled_dmr_counts":
        threshold = thresholds[0]
        source = mvi_root / "vmr" / ("var_%s/shared/input/selected_cells.allc.tsv" % threshold)
        bed = mvi_root / "vmr_dmr/shared/dmr/all_unique_hypo_DMRs.merged.bed"
        output = mvi_root / "vmr_dmr/shared/dmr/input.h5mu"
        features = sum(bool(line.strip()) for line in bed.open())
        run([methylvi_python, scripts / "Methylvi/vmr/02_build_vmr_methylvi_input.py",
             "--bed", bed, "--allc-table", source, "--annotation", annotation_table,
             "--work-dir", mvi_root / "vmr_dmr/shared/dmr/count_rows", "--output", output,
             "--threads", env["SCMO_THREADS"], "--min-covered-percent", "0", "--target-features", features,
             "--mc-context", allcools.get("mc_context", "CGN")])
        record([output])
    else:
        raise WorkflowError("no built-in task adapter for %s" % task)
    return 0


def vmr_environment(base, root, meth_root, annotation, chrom_sizes, blacklist,
                    methylvi_prefix, threshold, shared, maximum, max_features, expected, context):
    values = dict(base)
    values.update({
        "VMR_METHSCAN_RUN_DIR": str(meth_root), "VMR_METHSCAN_VARIANCE": str(threshold),
        "VMR_INPUT_MANIFEST": str(meth_root / "00_scanpy_selected/input_manifest.tsv"),
        "VMR_FILTERED_CELL_IDS": str(meth_root / "03_filtered/column_header.txt"),
        "VMR_SOURCE_BED": str(meth_root / "04_scan" / ("var_%s" % threshold) / "VMRs.bed"),
        "VMR_EXISTING_ALLC_DIR": str(shared / "input/allc"), "VMR_CHROM_SIZES": str(chrom_sizes),
        "VMR_BLACKLIST": str(blacklist), "VMR_ANNOTATION": str(annotation),
        "VMR_RESULTS_ROOT": str(shared), "VMR_INPUT_DIR": str(shared / "input"),
        "VMR_FILTERED_BED": str(shared / "input/vmr_filtered.bed"),
        "VMR_ALLC_TABLE": str(shared / "input/selected_cells.allc.tsv"),
        "VMR_MVI_INPUT": str(maximum / "input.h5mu"), "VMR_COUNT_ROWS": str(shared / "count_rows"),
        "VMR_MVI_RESULTS": str(maximum / "results"), "VMR_METHYLVI_ENV": str(methylvi_prefix),
        "VMR_EXPECTED_CELLS": str(expected), "VMR_TARGET_FEATURES": str(max_features),
        "VMR_MC_CONTEXT": str(context),
    })
    return values


def train_model(root, scripts, python, route, base_env, runner):
    input_path, results = route / "input.h5mu", route / "results"
    if not input_path.is_file():
        raise RuntimeError("MethylVI input is absent: %s" % input_path)
    extra = {"SCMO_MVI_INPUT": str(input_path), "SCMO_MVI_RESULTS": str(results),
             "VMR_MVI_RESULTS": str(results)}
    runner([python, scripts / "Methylvi/shared/04_train_methylvi.py"], extra)
    runner([python, scripts / "Methylvi/vmr/03_plot_vmr_methylvi_umap.py",
            "--input", results / "methylvi_embedding.h5ad", "--output-dir", results], extra)
    runner([python, scripts / "Methylvi/shared/05_plot_supervised_umap.py"], extra)
    runner([python, scripts / "Methylvi/shared/07_plot_methylation_qc.py",
            "--input", input_path, "--output-dir", results / "qc"], extra)
    (route / "model.COMPLETE").touch()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
