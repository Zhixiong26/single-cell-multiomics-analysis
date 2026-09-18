#!/usr/bin/env python3
"""Shared configuration, manifest, hashing, and JSON helpers."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


SCHEMA_VERSION = 2
SUPPORTED_SCHEMA_VERSIONS = {1, 2}
SAMPLE_COLUMNS = (
    "sample_id", "condition", "batch", "include", "rna_path", "rna_format",
    "allc_root", "allc_glob", "cell_id_prefix", "allc_cell_id_regex",
    "allc_cell_id_replacement",
)
LEGACY_SAMPLE_COLUMNS = SAMPLE_COLUMNS[:9]
RNA_FORMATS = {"", "10x_mtx", "10x_zip", "10x_h5"}


class WorkflowError(RuntimeError):
    pass


def load_structured(path: Path) -> Dict[str, Any]:
    path = Path(path).resolve()
    text = path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore
    except ImportError:
        try:
            value = json.loads(text)
        except ValueError as exc:
            raise WorkflowError(
                "%s requires PyYAML unless it is JSON-compatible YAML: %s" % (path, exc)
            )
    else:
        value = yaml.safe_load(text)
    if not isinstance(value, dict):
        raise WorkflowError("configuration root must be a mapping: %s" % path)
    return value


def write_json(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def signature(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def resolve_path(root: Path, value: Any) -> Optional[Path]:
    if value in (None, ""):
        return None
    path = Path(os.path.expandvars(os.path.expanduser(str(value))))
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def project_files(project: Path) -> Dict[str, Path]:
    root = Path(project).resolve()
    return {
        "root": root,
        "project": root / "config" / "project.yaml",
        "analysis": root / "config" / "analysis.yaml",
        "samples": root / "config" / "samples.tsv",
        "environments": root / "config" / "environments.tsv",
        "scheduler": root / "config" / "scheduler.yaml",
    }


def load_project(project: Path) -> Dict[str, Any]:
    files = project_files(project)
    missing = [str(path) for key, path in files.items() if key != "root" and not path.is_file()]
    if missing:
        raise WorkflowError("missing project configuration: %s" % ", ".join(missing))
    cfg = load_structured(files["project"])
    analysis = load_structured(files["analysis"])
    scheduler = load_structured(files["scheduler"])
    for label, value in (("project", cfg), ("analysis", analysis), ("scheduler", scheduler)):
        if value.get("schema_version") not in SUPPORTED_SCHEMA_VERSIONS:
            raise WorkflowError("%s schema_version must be one of %s" %
                                (label, sorted(SUPPORTED_SCHEMA_VERSIONS)))
    cfg["analysis"] = analysis.get("analysis", {})
    cfg["scheduler"] = scheduler
    cfg["project_root"] = str(files["root"])
    cfg["config_files"] = {key: str(value) for key, value in files.items() if key != "root"}
    return cfg


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: Iterable[str], rows: Iterable[Dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in writer.fieldnames})


def load_samples(project: Path, require_paths: bool = True) -> List[Dict[str, Any]]:
    files = project_files(project)
    with files["samples"].open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = tuple(reader.fieldnames or ())
        if columns not in {SAMPLE_COLUMNS, LEGACY_SAMPLE_COLUMNS}:
            raise WorkflowError("samples.tsv columns must be schema-v1 or schema-v2 columns: %s" % "\t".join(SAMPLE_COLUMNS))
        result: List[Dict[str, Any]] = []
        seen = set()
        errors = []
        for line_no, raw in enumerate(reader, 2):
            sample_id = raw["sample_id"].strip()
            if not re.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", sample_id):
                errors.append("line %d has invalid sample_id" % line_no)
            if sample_id in seen:
                errors.append("line %d duplicates sample_id %s" % (line_no, sample_id))
            seen.add(sample_id)
            include = raw["include"].strip().lower() in {"1", "true", "yes"}
            if raw["rna_format"].strip() not in RNA_FORMATS:
                errors.append("line %d has unsupported rna_format" % line_no)
            row: Dict[str, Any] = dict(raw)
            row["allc_cell_id_regex"] = row.get("allc_cell_id_regex") or ""
            row["allc_cell_id_replacement"] = row.get("allc_cell_id_replacement") or ""
            row["sample_id"] = sample_id
            row["include"] = include
            row["cell_id_prefix"] = raw["cell_id_prefix"].strip() or sample_id
            for field in ("rna_path", "allc_root"):
                resolved = resolve_path(files["root"], raw[field].strip())
                row[field] = str(resolved) if resolved else ""
                if include and require_paths and resolved and not resolved.exists():
                    errors.append("line %d %s does not exist: %s" % (line_no, field, resolved))
            if include and not row["rna_path"] and not row["allc_root"]:
                errors.append("line %d must provide RNA and/or ALLC input" % line_no)
            if row["rna_path"] and not raw["rna_format"].strip():
                errors.append("line %d requires rna_format" % line_no)
            result.append(row)
    if not any(row["include"] for row in result):
        errors.append("samples.tsv contains no included sample")
    if errors:
        raise WorkflowError("\n".join(errors))
    return result


def load_environments(project: Path) -> List[Dict[str, str]]:
    rows = load_tsv(project_files(project)["environments"])
    required = ("stage", "python", "executable", "version_command", "required")
    if rows and tuple(rows[0].keys()) != required:
        raise WorkflowError("environments.tsv columns must be exactly: %s" % "\t".join(required))
    return rows


def git_commit(root: Path) -> Optional[str]:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(root), stderr=subprocess.DEVNULL
        ).decode("utf-8").strip()
    except (OSError, subprocess.SubprocessError):
        return None


def parse_memory_mb(value: Any) -> int:
    text = str(value).strip().upper()
    match = re.match(r"^(\d+(?:\.\d+)?)\s*([KMGT]?)B?$", text)
    if not match:
        raise WorkflowError("invalid memory value: %s" % value)
    number = float(match.group(1))
    unit = match.group(2)
    factor = {"": 1, "K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}[unit]
    return int(number * factor)


def command_exists(command: str) -> bool:
    from shutil import which
    return which(command) is not None
