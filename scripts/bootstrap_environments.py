#!/usr/bin/env python3
"""Discover or create isolated Conda environments required by a generated project."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from _common import (
    RUNLOG_END, RUNLOG_START, STAGE_CONTEXT_END, STAGE_CONTEXT_START, WorkflowError, load_samples,
    marked_region, project_files, write_json,
)


PROFILE_ORDER = ("analysis_core", "methscan", "methylvi")
PROFILES: Dict[str, Dict[str, Any]] = {
    "analysis_core": {
        "stages": ("orchestrator", "scanpy_allcools"),
        "spec": "analysis-core.yaml",
        "executable": "bin/python",
        "check": [
            "bin/python", "-c",
            "import ALLCools,anndata,harmonypy,leidenalg,scanpy,scrublet,yaml; "
            "print('analysis_core', scanpy.__version__, anndata.__version__)",
        ],
    },
    "methscan": {
        "stages": ("methscan",),
        "spec": "methscan.yaml",
        "executable": "bin/methscan",
        "check": ["bin/methscan", "--version"],
    },
    "methylvi": {
        "stages": ("methylvi",),
        "spec": "methylvi.yaml",
        "executable": "bin/python",
        "check": [
            "bin/python", "-c",
            "import anndata,mudata,scanpy,scvi,torch; from scvi.external import METHYLVI; "
            "print('methylvi', scvi.__version__, torch.__version__, scanpy.__version__)",
        ],
    },
}
ENV_COLUMNS = ("stage", "python", "executable", "version_command", "required")


def required_profiles(project: Path) -> List[str]:
    active = [row for row in load_samples(project, require_paths=False) if row["include"]]
    has_rna = any(row["rna_path"] for row in active)
    has_allc = any(row["allc_root"] for row in active)
    result = ["analysis_core"] if has_rna or has_allc else []
    if has_allc:
        result.extend(["methscan", "methylvi"])
    return result


def read_environment_rows(project: Path) -> List[Dict[str, str]]:
    path = project_files(project)["environments"]
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ENV_COLUMNS:
            raise WorkflowError("environments.tsv has an invalid header")
        return list(reader)


def prefix_from_row(row: Dict[str, str]) -> Optional[Path]:
    python = row.get("python", "").strip()
    if not python:
        return None
    path = Path(python).expanduser()
    return path.parent.parent.resolve() if path.parent.name == "bin" else None


def check_prefix(prefix: Path, profile: str, timeout: int = 180) -> Dict[str, Any]:
    definition = PROFILES[profile]
    command = [str(prefix / definition["check"][0])] + list(definition["check"][1:])
    executable = prefix / definition["executable"]
    result: Dict[str, Any] = {
        "profile": profile, "prefix": str(prefix), "command": command,
        "valid": False, "version": "",
    }
    if not executable.is_file() or not os.access(str(executable), os.X_OK):
        result["error"] = "required executable is absent: %s" % executable
        return result
    try:
        environment = os.environ.copy()
        environment["PYTHONNOUSERSITE"] = "1"
        output = subprocess.check_output(command, stderr=subprocess.STDOUT, timeout=timeout, env=environment)
    except (OSError, subprocess.SubprocessError) as exc:
        result["error"] = str(exc)
        return result
    result["valid"] = True
    result["version"] = output.decode("utf-8", "replace").strip().splitlines()[-1]
    return result


def find_manager(requested: str) -> str:
    names = [requested] if requested != "auto" else ["mamba", "conda"]
    for name in names:
        path = shutil.which(name)
        if path:
            return str(Path(path).resolve())
    raise WorkflowError("no Conda-compatible manager found; install conda or mamba, or pass --manager")


def discover_prefixes(manager: str, rows: Iterable[Dict[str, str]]) -> List[Path]:
    candidates = []
    for row in rows:
        prefix = prefix_from_row(row)
        if prefix:
            candidates.append(prefix)
    candidates.append(Path(sys.prefix).resolve())
    try:
        output = subprocess.check_output([manager, "info", "--json"], stderr=subprocess.STDOUT, timeout=60)
        info = json.loads(output.decode("utf-8"))
        candidates.extend(Path(value).resolve() for value in info.get("envs", []))
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    for base in (Path.home() / "miniconda3" / "envs", Path.home() / "miniforge3" / "envs"):
        if base.is_dir():
            candidates.extend(path.resolve() for path in base.iterdir() if path.is_dir())
    unique = []
    seen = set()
    for path in candidates:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def environment_rows_for(profile: str, prefix: Path, version: str) -> List[Dict[str, str]]:
    definition = PROFILES[profile]
    python = prefix / "bin" / "python"
    executable = prefix / definition["executable"]
    check = ["/usr/bin/env", "PYTHONNOUSERSITE=1", str(prefix / definition["check"][0])] + list(definition["check"][1:])
    return [{
        "stage": stage,
        "python": str(python),
        "executable": str(executable if stage != "orchestrator" else python),
        "version_command": shlex.join(check),
        "required": "1",
    } for stage in definition["stages"]]


def write_environment_rows(project: Path, rows: List[Dict[str, str]]) -> None:
    path = project_files(project)["environments"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ENV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in ENV_COLUMNS})


def write_environment_report(project: Path, result: Dict[str, Any]) -> None:
    lines = [
        "# Environment report / 环境报告",
    ]
    path = Path(project) / "Scripts" / "Environment" / "Report.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    # This report is a snapshot of the current provisioning state, so it is
    # rewritten whole rather than accumulated. Two generated regions survive that
    # rewrite: the context block naming the project, and a run-log region if the
    # project has one. Neither is this tool's to invent or to discard.
    existing = ""
    if path.is_file():
        try:
            existing = path.read_text(encoding="utf-8")
        except OSError:
            existing = ""
    context = marked_region(existing, STAGE_CONTEXT_START, STAGE_CONTEXT_END)
    run_log = marked_region(existing, RUNLOG_START, RUNLOG_END)
    if context:
        lines.extend(["", context])
    lines.extend([
        "",
        "- Status / 状态: `%s`" % result["status"],
        "- Mode / 模式: `%s`" % result["mode"],
        "- Manager / 管理器: `%s`" % result["manager"],
        "- Timestamp / 时间: `%s`" % result["created_at"], "",
        "| Profile | Action | Prefix |", "|---|---|---|",
    ])
    for action in result["actions"]:
        lines.append("| `%s` | `%s` | `%s` |" % (
            action["profile"], action["action"], action["prefix"].replace("|", "\\|"),
        ))
    lines.extend([
        "", "Existing compatible environments were reused read-only; created environments are isolated.",
        "兼容的已有环境仅只读复用；新建环境均为隔离环境。", "",
    ])
    if run_log:
        lines.extend([run_log, ""])
    path.write_text("\n".join(lines).rstrip("\n") + "\n", encoding="utf-8")


def bootstrap(
    project: Path,
    execute: bool = False,
    manager_name: str = "auto",
    prefix_root: Optional[Path] = None,
    discover: bool = True,
) -> Dict[str, Any]:
    root = Path(project).resolve()
    files = project_files(root)
    rows = read_environment_rows(root)
    required = required_profiles(root)
    manager = find_manager(manager_name)
    target_root = (prefix_root or (root / ".environments")).resolve()
    specs = root / "environment-specs"
    candidates = discover_prefixes(manager, rows) if discover else []
    actions = []
    selected: Dict[str, Dict[str, Any]] = {}

    def fail(message: str) -> None:
        failed = {
            "schema_version": 1, "project": str(files["root"]),
            "created_at": datetime.now(timezone.utc).isoformat(),
            "mode": "execute" if execute else "plan", "manager": manager,
            "required_profiles": required, "actions": actions,
            "status": "failed", "error": message,
        }
        evidence = root / ".workflow" / "environment-bootstrap"
        write_json(evidence / "failure.json", failed)
        write_environment_report(root, failed)
        raise WorkflowError(message)

    for profile in required:
        definition = PROFILES[profile]
        declared = [prefix_from_row(row) for row in rows if row.get("stage") in definition["stages"]]
        search = [item for item in declared if item] + candidates
        checks = []
        for prefix in search:
            checked = check_prefix(prefix, profile)
            checks.append(checked)
            if checked["valid"]:
                selected[profile] = checked
                actions.append({"profile": profile, "action": "reuse", "prefix": str(prefix)})
                break
        if profile in selected:
            continue
        target = target_root / profile.replace("_", "-")
        action = {
            "profile": profile, "action": "create", "prefix": str(target),
            "spec": str(specs / definition["spec"]), "candidate_checks": checks,
        }
        actions.append(action)
        if not execute:
            continue
        spec = Path(action["spec"])
        if not spec.is_file():
            fail("environment spec is absent: %s" % spec)
        if target.exists() and (not target.is_dir() or any(target.iterdir())):
            fail("refusing non-empty invalid environment prefix: %s" % target)
        target.parent.mkdir(parents=True, exist_ok=True)
        command = [manager, "env", "create", "--prefix", str(target), "--file", str(spec)]
        action["command"] = command
        creation_environment = os.environ.copy()
        creation_environment["CONDA_ALWAYS_YES"] = "true"
        package_cache = target_root / ".package-cache"
        package_cache.mkdir(parents=True, exist_ok=True)
        creation_environment["CONDA_PKGS_DIRS"] = str(package_cache)
        try:
            subprocess.check_call(command, cwd=str(root), env=creation_environment)
        except (OSError, subprocess.SubprocessError) as exc:
            action["error"] = str(exc)
            fail("environment creation failed for %s: %s" % (profile, exc))
        checked = check_prefix(target, profile)
        action["verification"] = checked
        if not checked["valid"]:
            fail("created environment failed verification for %s: %s" % (profile, checked.get("error")))
        selected[profile] = checked

    if execute:
        managed_stages = {stage for profile in required for stage in PROFILES[profile]["stages"]}
        updated = [row for row in rows if row.get("stage") not in managed_stages]
        for profile in required:
            checked = selected.get(profile)
            if not checked:
                fail("environment remains unresolved: %s" % profile)
            updated.extend(environment_rows_for(profile, Path(checked["prefix"]), checked["version"]))
        write_environment_rows(root, updated)

    result = {
        "schema_version": 1,
        "project": str(files["root"]),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mode": "execute" if execute else "plan",
        "manager": manager,
        "required_profiles": required,
        "actions": actions,
        "status": "complete" if execute else "planned",
    }
    evidence = root / ".workflow" / "environment-bootstrap"
    write_json(evidence / ("result.json" if execute else "plan.json"), result)
    write_environment_report(root, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="create missing environments and update environments.tsv")
    parser.add_argument("--manager", default="auto", help="conda/mamba executable name or absolute path")
    parser.add_argument("--prefix-root", type=Path, help="root for newly created isolated environments")
    parser.add_argument("--no-discovery", action="store_true", help="skip scanning existing Conda prefixes")
    args = parser.parse_args()
    result = bootstrap(args.project, args.execute, args.manager, args.prefix_root, not args.no_discovery)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
