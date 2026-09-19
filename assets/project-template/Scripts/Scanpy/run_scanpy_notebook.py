#!/usr/bin/env python3
"""Execute the project Scanpy workflow notebook and record the evidence it produced.

This runner contains no analysis logic. It translates `config/` into the JSON parameter
sidecar that `Scripts/Scanpy/Notebooks/scanpy_workflow.ipynb` reads, executes that notebook
with the interpreter declared for the `scanpy_allcools` stage, and writes
`completion_summary.json` describing what the notebook produced. Every analysis decision,
threshold, figure and label lives in the notebook.

Sidecar keys the notebook reads are listed in `PARAMETERS_KEYS`; the harness refuses to run
when the configured values cannot produce a well-formed sidecar (missing inputs, an
unreviewed baseline mapping, a palette shorter than the sample list).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NOTEBOOK_RELATIVE = Path("Scripts") / "Scanpy" / "Notebooks" / "scanpy_workflow.ipynb"
EXECUTED_NAME = "scanpy_workflow_executed.ipynb"
PARAMETERS_NAME = "scanpy_parameters.json"
SUMMARY_NAME = "completion_summary.json"
ARTIFACTS_NAME = "notebook_artifacts.json"
KERNEL_NAME = "scanpy_project"
DEFAULT_PALETTE = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")
PROBE_PACKAGES = ("scanpy", "anndata", "harmonypy", "scrublet", "leidenalg", "numpy",
                  "pandas", "nbformat", "nbclient", "ipykernel", "jupyter_client")

# Defaults mirror the authoritative notebook, so a project configuration that predates a key
# still produces the same analysis the notebook performs when it is run by hand.
DEFAULTS: Dict[str, Any] = {
    "min_genes": 200,
    "max_genes": 6000,
    "min_counts": 500,
    "max_mt_percent": 5.0,
    "min_cells_per_gene": 3,
    "scrublet_min_cells": 100,
    "doublet_rate_per_1000": 0.004,
    "target_sum": 10000,
    "n_hvg": 2000,
    "hvg_batch_key": "cohort",
    "regress_covariates": False,
    "pca_n_comps": 50,
    "pca_svd_solver": "arpack",
    "n_pcs": 30,
    "n_neighbors": 15,
    "harmony_max_iter": 20,
    "harmony_sigma_value": 0.1,
    "umap_min_dist": 0.5,
    "umap_spread": 1.0,
    "leiden_resolution": 0.8,
    "candidate_min_score_margin": 0.20,
    "candidate_top_markers_per_cluster": 5,
    "seed": 0,
    "group_key": "group",
}
PARAMETERS_KEYS = tuple(sorted(set(DEFAULTS) | {
    "project_root", "formal_output_dir", "output_stem", "run_kind", "iteration_id",
    "overwrite_data_outputs", "analysis_confirmed", "analysis_signature",
    "annotation_profile_signature", "inputs", "sample_palette", "mitochondrial_prefixes",
    "ribosomal_prefixes", "dotplot_markers", "extra_candidate_panels", "cell_type_order",
    "cluster_to_cell_type", "epithelial_clusters", "epithelial_groups", "rare_clusters",
    "rare_groups",
}))
ITERATION_CHARSET = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")


class RunnerError(RuntimeError):
    pass


def nested(mapping: Any, *keys: str) -> Any:
    value = mapping
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def sanitize(text: str) -> str:
    cleaned = "".join(character if character in ITERATION_CHARSET else "_" for character in str(text))
    return cleaned or "run"


def load_context(root: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    sys.path.insert(0, str(root / "tools"))
    from _common import load_project, load_samples, WorkflowError

    try:
        cfg = load_project(root)
        samples = load_samples(root)
    except WorkflowError as exc:
        raise RunnerError(str(exc))
    return cfg, samples


def load_annotation_profile(root: Path, annotation: Dict[str, Any]) -> Tuple[Dict[str, str], Optional[str]]:
    from _common import load_structured, resolve_path, WorkflowError

    path = resolve_path(root, annotation.get("profile"))
    if path is None:
        return {}, None
    if not path.is_file():
        raise RunnerError("annotation profile does not exist: %s" % path)
    try:
        profile = load_structured(path)
    except WorkflowError as exc:
        raise RunnerError(str(exc))
    annotations = profile.get("annotations") or {}
    mapping: Dict[str, str] = {}
    for cluster, value in annotations.items():
        label = value.get("cell_type") if isinstance(value, dict) else value
        if not label:
            raise RunnerError("annotation profile has no cell_type for cluster %s" % cluster)
        mapping[str(cluster)] = str(label)
    signature = profile.get("analysis_signature")
    return mapping, (str(signature) if signature else None)


def build_parameters(root: Path, cfg: Dict[str, Any], samples: List[Dict[str, Any]],
                     output_dir: Path, run_id: str, iteration_id: Optional[str] = None) -> Dict[str, Any]:
    from _common import scanpy_analysis_signature

    analysis = nested(cfg, "analysis", "scanpy") or {}
    notebook_cfg = analysis.get("notebook") or {}
    project = cfg.get("project") or {}
    markers = cfg.get("markers") or {}
    annotation = cfg.get("annotation") or {}

    inputs: Dict[str, Dict[str, str]] = {}
    for row in samples:
        if not row["include"] or not row["rna_path"]:
            continue
        label = row["cell_id_prefix"]
        if label in inputs:
            raise RunnerError("two included samples share the cell_id_prefix %s; RNA cell IDs "
                              "would collide and could not be joined with ALLC cell IDs" % label)
        inputs[label] = {"path": row["rna_path"], "format": row["rna_format"]}
    if not inputs:
        raise RunnerError("no included RNA sample to analyse")

    palette = [str(value) for value in (analysis.get("sample_palette") or DEFAULT_PALETTE)]
    if len(palette) < len(inputs):
        raise RunnerError("analysis.scanpy.sample_palette provides %d colours for %d samples"
                          % (len(palette), len(inputs)))

    run_kind = str(notebook_cfg.get("run_kind") or "candidate")
    if run_kind not in {"baseline", "candidate"}:
        raise RunnerError("analysis.scanpy.notebook.run_kind must be baseline or candidate")
    if run_kind == "candidate":
        configured = iteration_id or notebook_cfg.get("iteration_id") or run_id
        resolved_iteration: Optional[str] = sanitize(configured)
        if resolved_iteration != str(configured):
            print("iteration id %r is not path-safe; using %r" % (str(configured), resolved_iteration),
                  file=sys.stderr)
    else:
        if iteration_id or notebook_cfg.get("iteration_id"):
            raise RunnerError("baseline runs must not set an iteration id")
        resolved_iteration = None

    cluster_to_cell_type, profile_signature = load_annotation_profile(root, annotation)
    confirmed = bool(notebook_cfg.get("analysis_confirmed", False))
    if run_kind == "baseline" and not confirmed:
        raise RunnerError("baseline runs require analysis.scanpy.notebook.analysis_confirmed: true "
                          "after the parameters and the reviewed annotation have been checked")
    if run_kind == "baseline" and not cluster_to_cell_type:
        raise RunnerError("baseline runs require a reviewed annotation profile under "
                          "config/annotation.yaml; write one with tools/record_annotation_review.py")

    panels = {str(label): [str(gene) for gene in genes]
              for label, genes in (markers.get("dotplot_markers") or {}).items()}
    extra = {str(label): [str(gene) for gene in genes]
             for label, genes in (markers.get("extra_candidate_panels") or {}).items()}
    order = [str(label) for label in (markers.get("cell_type_order") or [])]
    if panels:
        if order and set(order) != set(panels):
            raise RunnerError("markers.cell_type_order must list exactly the markers.dotplot_markers "
                              "labels; the annotation dotplot requires one row per panel group")
        if cluster_to_cell_type and set(cluster_to_cell_type.values()) != set(panels):
            raise RunnerError("every mapped cell type must have a markers.dotplot_markers group; "
                              "mapped=%s panel=%s" % (sorted(set(cluster_to_cell_type.values())),
                                                      sorted(panels)))
    else:
        print("no markers.dotplot_markers configured: the candidate run will emit data-driven "
              "cluster evidence for annotation review", file=sys.stderr)

    parameters: Dict[str, Any] = {
        "project_root": str(root),
        "formal_output_dir": str(output_dir),
        "output_stem": str(notebook_cfg.get("output_stem") or project.get("id") or "scanpy"),
        "run_kind": run_kind,
        "iteration_id": resolved_iteration,
        "overwrite_data_outputs": bool(notebook_cfg.get("overwrite_data_outputs", True)),
        "analysis_confirmed": confirmed,
        "inputs": inputs,
        "sample_palette": palette,
        "mitochondrial_prefixes": [str(value) for value in
                                   (project.get("mitochondrial_prefixes") or ["MT-", "mt-"])],
        "ribosomal_prefixes": [str(value) for value in
                               (project.get("ribosomal_prefixes") or ["RPL", "RPS"])],
        "group_key": str(nested(analysis, "group_key") or DEFAULTS["group_key"]),
        "dotplot_markers": panels,
        "extra_candidate_panels": extra,
        "cell_type_order": order,
        "cluster_to_cell_type": cluster_to_cell_type,
        "epithelial_clusters": [str(value) for value in (markers.get("epithelial_clusters") or [])],
        "epithelial_groups": [str(value) for value in (markers.get("epithelial_groups") or [])],
        "rare_clusters": [str(value) for value in (markers.get("rare_clusters") or [])],
        "rare_groups": [str(value) for value in (markers.get("rare_groups") or [])],
        "annotation_profile_signature": profile_signature,
    }
    qc_keys = {"min_genes": ("qc", "min_genes"), "max_genes": ("qc", "max_genes"),
               "min_counts": ("qc", "min_counts"), "max_mt_percent": ("qc", "max_mt_percent"),
               "min_cells_per_gene": ("qc", "min_cells_per_gene"),
               "scrublet_min_cells": ("scrublet", "min_cells"),
               "doublet_rate_per_1000": ("scrublet", "doublet_rate_per_1000"),
               "target_sum": ("target_sum",), "n_hvg": ("hvg", "n_top_genes"),
               "hvg_batch_key": ("hvg", "batch_key"),
               "regress_covariates": ("regress_covariates",),
               "pca_n_comps": ("pca", "n_comps"), "pca_svd_solver": ("pca", "svd_solver"),
               "n_pcs": ("neighbors", "n_pcs"), "n_neighbors": ("neighbors", "n_neighbors"),
               "harmony_max_iter": ("harmony", "max_iter"),
               "harmony_sigma_value": ("harmony", "sigma"),
               "umap_min_dist": ("umap", "min_dist"), "umap_spread": ("umap", "spread"),
               "leiden_resolution": ("leiden", "resolution"),
               "candidate_min_score_margin": ("candidate", "min_score_margin"),
               "candidate_top_markers_per_cluster": ("candidate", "top_markers_per_cluster"),
               "seed": ("seed",)}
    for key, path in qc_keys.items():
        value = nested(analysis, *path)
        parameters[key] = DEFAULTS[key] if value is None else value

    parameters["analysis_signature"] = scanpy_analysis_signature(parameters)
    if run_kind == "baseline" and profile_signature != parameters["analysis_signature"]:
        # Fail before the analysis runs: a profile recorded for other parameters must never be
        # reused, and the notebook refuses it too.
        raise RunnerError("the annotation profile was recorded for a different analysis signature "
                          "(%s versus %s); re-run tools/record_annotation_review.py for the current "
                          "parameters" % (profile_signature or "absent", parameters["analysis_signature"]))
    missing = [key for key in PARAMETERS_KEYS if key not in parameters]
    if missing:  # a sidecar without a notebook key would fail inside the kernel instead of here
        raise RunnerError("sidecar is missing keys the notebook reads: %s" % ", ".join(sorted(missing)))
    extra_keys = sorted(set(parameters) - set(PARAMETERS_KEYS))
    if extra_keys:
        raise RunnerError("sidecar carries keys the notebook does not read: %s" % ", ".join(extra_keys))
    return parameters


def kernel_environment(output_dir: Path, python: Path, sidecar: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Build the kernel spec and the process environment the kernel inherits.

    `JUPYTER_PATH` makes the kernel spec below the only one this run can select, so the kernel is
    always the `scanpy_allcools` interpreter; `SCMO_SCANPY_PARAMETERS` is how the notebook learns
    the sidecar it must read.
    """
    kernel_root = output_dir / ".jupyter"
    spec_dir = kernel_root / "kernels" / KERNEL_NAME
    spec_dir.mkdir(parents=True, exist_ok=True)
    spec = {"argv": [str(python), "-m", "ipykernel_launcher", "-f", "{connection_file}"],
            "display_name": "Scanpy project environment", "language": "python",
            "metadata": {"debugger": False}}
    (spec_dir / "kernel.json").write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n",
                                          encoding="utf-8")
    environment = os.environ.copy()
    environment["JUPYTER_PATH"] = str(kernel_root)
    environment["SCMO_SCANPY_PARAMETERS"] = str(sidecar)
    return environment, spec


def probe_versions(python: Path) -> Dict[str, str]:
    code = "\n".join([
        "import importlib, json",
        "names = %r" % (PROBE_PACKAGES,),
        "versions = {}",
        "for name in names:",
        "    try:",
        "        versions[name] = str(getattr(importlib.import_module(name), '__version__', 'unknown'))",
        "    except Exception as exc:",
        "        versions[name] = 'unavailable (%s)' % type(exc).__name__",
        "print(json.dumps(versions, sort_keys=True))",
    ])
    try:
        output = subprocess.check_output([str(python), "-c", code], stderr=subprocess.DEVNULL, timeout=300)
    except (OSError, subprocess.SubprocessError) as exc:
        return {"probe": "failed (%s)" % type(exc).__name__}
    for line in reversed(output.decode("utf-8", "replace").splitlines()):
        try:
            return json.loads(line)
        except ValueError:
            continue
    return {"probe": "unreadable"}


def execute(root: Path, output_dir: Path, notebook_path: Path, environment: Dict[str, str],
            timeout: Optional[float]) -> Tuple[str, Optional[str]]:
    try:
        import nbformat
        from nbclient import NotebookClient
        from nbclient.exceptions import CellExecutionError, CellTimeoutError, DeadKernelError
    except ImportError as exc:
        raise RunnerError("the scanpy_allcools interpreter lacks the notebook runtime (%s); add "
                          "ipykernel, nbclient and nbformat to environment-specs/analysis-core.yaml "
                          "and re-run tools/bootstrap_environments.py" % exc)

    os.environ.update(environment)
    notebook = nbformat.read(str(notebook_path), as_version=4)
    client = NotebookClient(notebook, timeout=timeout, kernel_name=KERNEL_NAME,
                            resources={"metadata": {"path": str(root)}},
                            allow_errors=False, record_timing=True)

    def announce(**kwargs: Any) -> None:
        cell = kwargs.get("cell") or {}
        if cell.get("cell_type") != "code":
            return
        source = "".join(cell.get("source", ""))
        first = next((value.strip() for value in source.splitlines()
                      if value.strip() and not value.strip().startswith("#")), "")
        print("[cell %s] %s" % (kwargs.get("cell_index"), first[:110]), flush=True)

    def forward(**kwargs: Any) -> None:
        for output in (kwargs.get("cell") or {}).get("outputs", []):
            text = "".join(output.get("text", "")) if output.get("output_type") == "stream" else ""
            if text:
                sys.stderr.write(text)
        sys.stderr.flush()

    client.on_cell_start = announce
    client.on_cell_complete = forward

    failure: Optional[str] = None
    try:
        client.execute()
        status = "succeeded"
    except CellExecutionError as exc:  # a cell raised: its traceback is the evidence
        status = "failed"
        failure = str(exc)[-6000:]
        print("notebook cell failed:\n%s" % failure, file=sys.stderr)
    except (CellTimeoutError, DeadKernelError) as exc:  # the run cannot be trusted past this point
        status = "failed"
        failure = "%s: %s" % (type(exc).__name__, exc)
        print("notebook execution aborted: %s" % failure, file=sys.stderr)
    except Exception as exc:  # kernel startup, unreadable notebook, unwritable outputs
        status = "failed"
        failure = "%s: %s" % (type(exc).__name__, exc)
        print("notebook execution failed: %s" % failure, file=sys.stderr)
    finally:
        executed_path = output_dir / EXECUTED_NAME
        nbformat.write(client.nb, str(executed_path))
    return status, failure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--iteration-id", default=None)
    parser.add_argument("--timeout-seconds", type=float, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg, samples = load_context(root)

    template = root / NOTEBOOK_RELATIVE
    if not template.is_file():
        raise RunnerError("workflow notebook is absent: %s" % template)

    from _common import load_environments, resolve_path, sha256_file, write_json

    environment_row = next((row for row in load_environments(root)
                            if row["stage"].strip() == "scanpy_allcools"), None)
    if environment_row is None:
        raise RunnerError("config/environments.tsv declares no scanpy_allcools stage")
    python = resolve_path(root, environment_row["python"])
    if python is None or not python.is_file():
        raise RunnerError("the scanpy_allcools Python is absent: %s" % environment_row["python"])

    notebook_cfg = nested(cfg, "analysis", "scanpy", "notebook") or {}
    timeout = args.timeout_seconds if args.timeout_seconds is not None else notebook_cfg.get("timeout_seconds")
    parameters = build_parameters(root, cfg, samples, output_dir, args.run_id, args.iteration_id)

    sidecar = output_dir / PARAMETERS_NAME
    write_json(sidecar, parameters)
    kernel_env, spec = kernel_environment(output_dir, python, sidecar)
    started = datetime.now(timezone.utc)
    clock = time.time()
    status, failure = execute(root, output_dir, template, kernel_env, timeout)
    duration = round(time.time() - clock, 3)

    artifacts: List[str] = []
    reported: Dict[str, Any] = {}
    search = [output_dir]
    if parameters["iteration_id"]:
        search.insert(0, output_dir / "iterations" / str(parameters["iteration_id"]))
    artifacts_path = next((folder / ARTIFACTS_NAME for folder in search
                           if (folder / ARTIFACTS_NAME).is_file()), None)
    if artifacts_path is None:
        artifacts = []
        if status == "succeeded":
            failure = "the notebook completed but wrote no %s" % ARTIFACTS_NAME
            status = "failed"
    else:
        # The notebook reports absolute paths it verified on disk; the harness only checks them.
        reported = json.loads(artifacts_path.read_text(encoding="utf-8"))
        artifacts = [str(Path(value).resolve()) for value in reported.get("artifacts", [])]

    fixed = [output_dir / EXECUTED_NAME, sidecar]
    summary = {
        "task": "scanpy", "status": status, "run_kind": parameters["run_kind"],
        "iteration_id": parameters["iteration_id"], "run_id": args.run_id,
        "started_at": started.isoformat(), "duration_seconds": duration,
        "project": str(root), "output_dir": str(output_dir),
        "notebook": str((output_dir / EXECUTED_NAME).resolve()),
        "notebook_template": {"path": str(template.resolve()), "sha256": sha256_file(template)},
        "parameters": str(sidecar.resolve()),
        "analysis_signature": parameters["analysis_signature"],
        "annotation_profile_signature": parameters["annotation_profile_signature"],
        "kernel": {"name": KERNEL_NAME, "python": str(python), "argv": spec["argv"],
                   "environment": environment_row.get("version_command", "")},
        "packages": probe_versions(python),
        "annotation_status": reported.get("annotation_status"),
        "annotation_key": reported.get("annotation_column"),
        "annotation_table": reported.get("annotation_table"),
        "n_cells": reported.get("n_cells"), "n_genes": reported.get("n_genes"),
        "n_clusters": reported.get("n_clusters"),
        "n_figures": reported.get("n_figures"),
        "artifacts": sorted({str(Path(value).resolve()) for value in artifacts + [str(path) for path in fixed]}),
    }
    if failure:
        summary["failure"] = failure
    write_json(output_dir / SUMMARY_NAME, summary)

    print(json.dumps({"status": status, "output_dir": str(output_dir),
                      "executed_notebook": summary["notebook"],
                      "annotation_status": summary["annotation_status"],
                      "n_cells": summary["n_cells"], "n_clusters": summary["n_clusters"],
                      "n_figures": reported.get("n_figures")},
                     indent=2, sort_keys=True), flush=True)
    if status != "succeeded":
        print("scanpy notebook failed; inspect %s" % summary["notebook"], file=sys.stderr)
        return 1
    if parameters["run_kind"] == "candidate":
        print("candidate run finished; review the cluster evidence, then record the annotation and "
              "re-run with run_kind: baseline:\n"
              "  tools/record_annotation_review.py --project . --run-id %s --worksheet review.tsv\n"
              "  tools/record_annotation_review.py --project . --run-id %s --mapping review.tsv"
              % (args.run_id, args.run_id), file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RunnerError, OSError, ValueError, KeyError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
