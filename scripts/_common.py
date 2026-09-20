#!/usr/bin/env python3
"""Shared configuration, manifest, hashing, and JSON helpers."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


SCHEMA_VERSION = 2
SUPPORTED_SCHEMA_VERSIONS = {1, 2}
SAMPLE_COLUMNS = (
    "sample_id", "condition", "batch", "include", "rna_path", "rna_format",
    "allc_root", "allc_glob", "cell_id_prefix", "allc_cell_id_regex",
    "allc_cell_id_replacement",
)
LEGACY_SAMPLE_COLUMNS = SAMPLE_COLUMNS[:9]
RNA_FORMATS = {"", "10x_mtx", "10x_zip", "10x_h5"}
BUILTIN_TASKS = {
    "scanpy", "methscan_select_convert", "methscan_prepare", "methscan_filter",
    "methscan_smooth", "methscan_pairwise_dmr", "methscan_hypo_heatmaps",
    "methscan_pooled_dmr", "allcools_features", "pooled_dmr_prepare",
    "pooled_dmr_counts", "workflow_summary",
}
BUILTIN_TASK_PREFIXES = (
    "methscan_vmr_", "methylvi_allcools_", "methylvi_vmr_features_",
    "methylvi_vmr_dmr_", "methylvi_vmr_",
)

# Scheduler and wrapper states that mean a task will never finish on its own.
TERMINAL_BAD_STATES = {
    "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE", "REVOKED", "SPECIAL_EXIT", "failed",
}

# The root Report.md run log. One region holds every per-run record; inside it
# each run gets its own marker pair so a re-inspection updates that run in place
# instead of appending a duplicate. Only the region is tool-owned: the heading
# above it and anything outside it is human text and is preserved verbatim.
RUNLOG_START = "<!-- SCMO-RUNLOG:START -->"
RUNLOG_END = "<!-- SCMO-RUNLOG:END -->"
RUNLOG_TITLE = "## Run log / 运行记录"
RUNLOG_FORMAT = 1
# Superseded single-block format written by older releases; migrated on first rewrite.
LEGACY_RUN_STATE_START = "<!-- SCMO-RUN-STATE:START -->"
LEGACY_RUN_STATE_END = "<!-- SCMO-RUN-STATE:END -->"

# Stage reports. Each stage's Report.md accumulates the same run summaries as the
# root Report, filtered to the tasks that stage owns, so a reader of the Scanpy
# report does not have to read the MethylVI runs around it. The root Report stays
# the complete record; these are per-stage views of one corpus of summaries, not
# a second source of truth. Only reports a project already has are rewritten --
# the file's human half is the stage's contract and cannot be invented here.
STAGE_REPORTS = (
    ("Scanpy", "Scripts/Scanpy/Report.md"),
    ("Methscan", "Scripts/Methscan/Report.md"),
    ("Methylvi/allcools", "Scripts/Methylvi/allcools/Report.md"),
    ("Methylvi", "Scripts/Methylvi/Report.md"),
    ("Methylvi/vmr_dmr", "Scripts/Methylvi/vmr_dmr/Report.md"),
)

# Which task belongs to which stage report. Prefix order matters: the
# `methylvi_vmr_dmr_` trainers would otherwise match the broader `methylvi_vmr_`.
STAGE_TASKS = {
    "scanpy": "Scanpy",
    "methscan_select_convert": "Methscan",
    "methscan_prepare": "Methscan",
    "methscan_filter": "Methscan",
    "methscan_smooth": "Methscan",
    "methscan_pairwise_dmr": "Methscan",
    "methscan_hypo_heatmaps": "Methscan",
    "methscan_pooled_dmr": "Methscan",
    "allcools_features": "Methylvi/allcools",
    "pooled_dmr_prepare": "Methylvi/vmr_dmr",
    "pooled_dmr_counts": "Methylvi/vmr_dmr",
}
STAGE_TASK_PREFIXES = (
    ("methscan_vmr_", "Methscan"),
    ("methylvi_vmr_dmr_", "Methylvi/vmr_dmr"),
    ("methylvi_vmr_", "Methylvi"),
    ("methylvi_allcools_", "Methylvi/allcools"),
)

# The generated context block in every stage document. It names the project the
# copy belongs to, so the reference project's example sample names and worked
# numbers in the template prose are never mistaken for this project's own. It is
# written once by init_project.py and preserved by every tool that rewrites part
# of a stage document, including bootstrap_environments.py, whose report is a
# whole-file snapshot rather than a run log.
STAGE_CONTEXT_START = "<!-- SCMO-STAGE-CONTEXT:START -->"
STAGE_CONTEXT_END = "<!-- SCMO-STAGE-CONTEXT:END -->"


def marked_region(text: str, start: str, end: str) -> Optional[str]:
    """The block from `start` through `end` inclusive, or None when either is absent.

    Used to carry a region across a whole-file rewrite. Markers are matched as
    they appear in the file rather than line-anchored, because these blocks are
    written by a single producer and a partial match is a corrupt file that the
    caller should leave alone rather than half-preserve.
    """
    begin = text.find(start)
    if begin < 0:
        return None
    finish = text.find(end, begin + len(start))
    if finish < 0:
        return None
    return text[begin:finish + len(end)]


def stage_for_task(task_id: Any) -> Optional[str]:
    """The stage report a task belongs to, or None when no stage owns it.

    `workflow_summary` is deliberately unowned: it summarises the whole run and
    would otherwise appear in every stage log as though each stage had run it.
    """
    text = str(task_id or "")
    if text in STAGE_TASKS:
        return STAGE_TASKS[text]
    for prefix, stage in STAGE_TASK_PREFIXES:
        if text.startswith(prefix):
            return stage
    return None


def _stage_status(tasks: Sequence[Dict[str, Any]]) -> str:
    """A stage's own verdict, from its tasks rather than from the run's status.

    A run can complete every stage but one; reporting the run's status here would
    let that stage's report read `complete` while its own task failed.
    """
    if any(str(row.get("state") or "") in TERMINAL_BAD_STATES for row in tasks):
        return "failed"
    if all(row.get("evidence_valid") for row in tasks):
        return "complete"
    return "unfinished"

RUN_MARKER_RE = re.compile(r"^[ \t]*<!--[ \t]*SCMO-RUN:([^\s<>]+):(START|END)[ \t]*-->[ \t]*$")
# Recognises only a whole line that is the region marker itself, so a run record
# that merely mentions the marker in prose (or in backticks) is left alone.
RUNLOG_MARKER_RE = re.compile(r"^[ \t]*<!--[ \t]*SCMO-RUNLOG:(?:START|END)[ \t]*-->[ \t]*$")
# A marker key may not contain "<" or ">" or whitespace, so it can never break
# out of the HTML comment. run_id values are not that constrained (plan_workflow
# only rejects "", "/", ".", ".."), so unsafe ids are replaced by a digest.
RUN_KEY_SAFE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
RUN_HEADING_RE = re.compile(r"^###[ \t]+`([^`]+)`")
RUN_CHECKED_RE = re.compile(r"^[ \t]*-[ \t]*Checked[^\n`]*`([^`]+)`")


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


def write_text_atomic(path: Path, text: str) -> None:
    # Mirror write_json: a crash mid-write must not truncate a document that is
    # the only copy of a project's run history. The dot prefix keeps a crashed
    # temporary file out of sight in the project root.
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name("." + path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def as_bool(value: Any) -> bool:
    """The one truthiness rule for a declared include flag.

    The project README renders what a sample's `include` means, so it must use
    the same rule `load_samples` uses; otherwise the document can claim a sample
    is excluded while the pipeline includes it.
    """
    return str(value).strip().lower() in {"1", "true", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def signature(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


SCANPY_PLUMBING_KEYS = frozenset({
    "project_root", "formal_output_dir", "output_stem", "run_kind", "iteration_id",
    "overwrite_data_outputs", "analysis_confirmed", "analysis_signature",
    "annotation_profile_signature", "sample_palette", "group_key",
    "dotplot_markers", "extra_candidate_panels", "cell_type_order", "cluster_to_cell_type",
    "epithelial_clusters", "epithelial_groups", "rare_clusters", "rare_groups",
    "candidate_min_score_margin", "candidate_top_markers_per_cluster",
})
"""Sidecar keys excluded from the Scanpy analysis signature.

Everything else in the sidecar is treated as analysis-relevant, so a parameter added to the
notebook invalidates a recorded annotation profile by default instead of silently escaping
the signature. The excluded keys only steer paths, run identity, figure styling, and the
reviewed label set; `candidate_*` keys shape the candidate audit table, which is explicitly
non-authoritative and never feeds the reviewed cluster mapping.
"""


def notebook_cells_digest(path: Path) -> str:
    """Canonical digest of a notebook's cell types and sources, ignoring outputs and metadata.

    Executed notebooks carry outputs, execution counts and timings that differ on every run;
    the digest deliberately covers only the inputs that determine the analysis.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    cells = payload.get("cells", [])
    if not isinstance(cells, list) or not cells:
        raise WorkflowError("notebook has no cells: %s" % path)
    canonical = [{"cell_type": cell.get("cell_type", ""), "source": "".join(cell.get("source", []))}
                 for cell in cells]
    return signature(canonical)


def scanpy_analysis_signature(parameters: Dict[str, Any]) -> str:
    """Signature of the analysis-relevant subset of a Scanpy run parameter sidecar."""
    if not isinstance(parameters, dict):
        raise WorkflowError("scanpy parameters must be a mapping")
    relevant = {key: value for key, value in parameters.items() if key not in SCANPY_PLUMBING_KEYS}
    return signature(relevant)


def task_override(commands: Dict[str, Any], task_id: str) -> Any:
    if task_id in commands:
        return commands[task_id]
    for pattern, command in commands.items():
        if pattern.endswith("*") and task_id.startswith(pattern[:-1]):
            return command
    return None


def task_is_implemented(task_id: str, commands: Dict[str, Any]) -> bool:
    return (task_override(commands, task_id) is not None
            or task_id in BUILTIN_TASKS
            or any(task_id.startswith(prefix) for prefix in BUILTIN_TASK_PREFIXES))


def validate_recorded_outputs(values: Iterable[Any], allow_empty: bool = False) -> tuple[bool, str | None]:
    outputs = [Path(str(value)) for value in values]
    if not outputs and not allow_empty:
        return False, "task declares no output evidence"
    for path in outputs:
        if not path.exists():
            return False, "declared task output is absent: %s" % path
        if path.name == "task_outputs.json":
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                artifacts = payload.get("artifacts", [])
            except (OSError, ValueError, TypeError, AttributeError):
                return False, "task_outputs.json is not a valid JSON object: %s" % path
            if not isinstance(artifacts, list) or not artifacts:
                return False, "task_outputs.json must contain a non-empty artifacts list: %s" % path
            relative = [str(artifact) for artifact in artifacts if not Path(str(artifact)).is_absolute()]
            if relative:
                return False, "recorded artifact paths must be absolute: %s" % ", ".join(relative)
            missing = [str(artifact) for artifact in artifacts if not Path(str(artifact)).exists()]
            if missing:
                return False, "recorded artifacts are absent: %s" % ", ".join(missing)
    return True, None


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


# Input data lives under the project's own Data/ directory. A declared source is
# either a path to link or a URL to fetch; both end up as Data/<name>, so a reader
# can see every input a project consumes without reading the sample manifest.
DATA_ROOT_NAME = "Data"
DATA_URL_SCHEMES = ("https", "http", "ftp")
DATA_STATE_ORDER = ("ready", "link", "dangling", "present", "download", "conflict")


def data_root(root: Path) -> Path:
    return Path(root) / DATA_ROOT_NAME


def data_sources(block: Any, base: Path) -> List[Dict[str, Any]]:
    """Normalize and validate a `data.sources` block into linkable entries.

    `base` anchors a relative `path` so one intake resolves to the same targets
    no matter which directory it is run from.

    A `url` source must carry a sha256. Every later stage reads whatever sits at
    the destination as data, and a truncated transfer is indistinguishable from a
    complete one by size alone, so an unchecked download is not accepted.
    """
    if block in (None, ""):
        return []
    if not isinstance(block, dict):
        raise WorkflowError("data must be a mapping with a sources list")
    raw_sources = block.get("sources") or []
    if not isinstance(raw_sources, list):
        raise WorkflowError("data.sources must be a list")
    entries: List[Dict[str, Any]] = []
    errors: List[str] = []
    seen = set()
    for index, raw in enumerate(raw_sources, 1):
        if not isinstance(raw, dict):
            errors.append("data.sources[%d] must be a mapping" % index)
            continue
        name = str(raw.get("name") or "").strip()
        if not name or name in {".", ".."} or "/" in name or name.startswith("."):
            errors.append("data.sources[%d] has an invalid name: it must be a single path component "
                          "without a leading dot, not %r" % (index, raw.get("name")))
            continue
        if name in seen:
            errors.append("data.sources[%d] duplicates the name %s" % (index, name))
            continue
        seen.add(name)
        path = str(raw.get("path") or "").strip()
        url = str(raw.get("url") or "").strip()
        if bool(path) == bool(url):
            errors.append("data.sources[%s] must declare exactly one of path or url" % name)
            continue
        entry: Dict[str, Any] = {"name": name, "target": "", "url": "", "sha256": ""}
        if path:
            entry["target"] = str(resolve_path(base, path))
        else:
            scheme = url.split("://", 1)[0].lower() if "://" in url else ""
            if scheme not in DATA_URL_SCHEMES:
                errors.append("data.sources[%s] url must start with one of %s: %s"
                              % (name, ", ".join(DATA_URL_SCHEMES), url))
                continue
            entry["url"] = url
            entry["sha256"] = str(raw.get("sha256") or "").strip().lower()
            if not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]):
                errors.append("data.sources[%s] is fetched from a url and therefore requires a "
                              "64-hex sha256 of the finished file" % name)
                continue
        entries.append(entry)
    if errors:
        raise WorkflowError("\n".join(errors))
    return entries


def stored_data_sources(cfg: Any) -> List[Dict[str, Any]]:
    """Read the normalized `data.sources` a generated project stores.

    Validated rather than trusted: editing `config/project.yaml` by hand is the
    normal way to retarget a project, so it gets the same treatment an intake does
    and the caller sees one entry shape either way.
    """
    block = (cfg or {}).get("data") if isinstance(cfg, dict) else None
    if not block:
        return []
    if not isinstance(block, dict) or not isinstance(block.get("sources"), list):
        raise WorkflowError("data.sources must be a list")
    entries: List[Dict[str, Any]] = []
    errors: List[str] = []
    seen = set()
    for index, raw in enumerate(block["sources"], 1):
        if not isinstance(raw, dict):
            errors.append("data.sources[%d] must be a mapping" % index)
            continue
        entry = {"name": str(raw.get("name") or "").strip(), "target": str(raw.get("target") or ""),
                 "url": str(raw.get("url") or ""), "sha256": str(raw.get("sha256") or "").strip().lower()}
        if not entry["name"] or entry["name"] in {".", ".."} or "/" in entry["name"] or entry["name"].startswith("."):
            errors.append("data.sources[%d] has an invalid name: %r" % (index, raw.get("name")))
            continue
        if entry["name"] in seen:
            errors.append("data.sources[%d] duplicates the name %s" % (index, entry["name"]))
            continue
        seen.add(entry["name"])
        if bool(entry["target"]) == bool(entry["url"]):
            errors.append("data.sources[%s] must declare exactly one of target or url" % entry["name"])
            continue
        if entry["target"] and not Path(entry["target"]).is_absolute():
            errors.append("data.sources[%s] target must be an absolute path" % entry["name"])
            continue
        if entry["url"] and not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]):
            errors.append("data.sources[%s] is fetched from a url and therefore requires a 64-hex sha256"
                          % entry["name"])
            continue
        entries.append(entry)
    if errors:
        raise WorkflowError("\n".join(errors))
    return entries


def plan_data_sources(root: Path, entries: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Say what each declared source needs, without touching the filesystem.

    `dangling` means the entry cannot be read from here, whichever way it got that
    way: a link that is absent with a missing target, or a link that is already
    exactly right but points at a path this host cannot see. Keeping that distinct
    from `ready` is what stops a correctly-shaped link to a vanished mount from
    being reported as nothing to do.

    `conflict` is deliberately its own state rather than an error at plan time:
    something that is not the declared entry is already in place, and deciding
    which copy of the data is real is a human call, not one this plan can make.
    """
    plans: List[Dict[str, Any]] = []
    for entry in entries:
        destination = data_root(root) / entry["name"]
        target_exists = Path(entry["target"]).exists() if entry["target"] else False
        plan = dict(entry, destination=str(destination), target_exists=target_exists)
        if destination.is_symlink():
            if not (entry["target"] and destination.resolve() == Path(entry["target"])):
                plan["state"] = "conflict"
            else:
                plan["state"] = "ready" if target_exists else "dangling"
        elif destination.exists():
            plan["state"] = "present" if entry["url"] else "conflict"
        elif entry["target"]:
            plan["state"] = "link" if target_exists else "dangling"
        else:
            plan["state"] = "download"
        plans.append(plan)
    return plans


def apply_data_sources(root: Path, plans: Sequence[Dict[str, Any]], download: bool = False) -> List[Dict[str, Any]]:
    """Carry out a plan, refusing anything the plan could not resolve itself.

    An unresolved link is still written: the target may be reachable only from the
    execution context, and validation there names it precisely. Refusing instead
    would make it impossible to generate a project on a host that cannot see the
    compute filesystem.
    """
    results: List[Dict[str, Any]] = []
    for plan in plans:
        destination = Path(plan["destination"])
        state = plan["state"]
        if state == "ready":
            results.append(dict(plan, outcome="unchanged"))
            continue
        if state == "conflict":
            raise WorkflowError("%s already exists and is not the entry this source declares; "
                                "move it aside before linking %s" % (destination, plan["name"]))
        if state == "dangling" and destination.is_symlink():
            # The link is already what was asked for; only its target is out of
            # reach here, and re-creating it would fail with EEXIST for no gain.
            results.append(dict(plan, outcome="unresolved"))
            continue
        if state == "download" and not download:
            results.append(dict(plan, outcome="pending"))
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        if state in {"link", "dangling"}:
            destination.symlink_to(plan["target"])
            results.append(dict(plan, outcome="linked" if state == "link" else "linked-unresolved"))
        elif state == "download":
            observed = download_source(destination, plan["url"], plan["sha256"])
            results.append(dict(plan, outcome="downloaded", observed_sha256=observed))
        else:
            observed = sha256_file(destination)
            if observed.lower() != plan["sha256"]:
                raise WorkflowError("%s exists but its sha256 is %s, not the declared %s; "
                                    "move it aside and fetch it again" % (destination, observed, plan["sha256"]))
            results.append(dict(plan, outcome="verified", observed_sha256=observed))
    return results


def download_source(destination: Path, url: str, expected_sha256: str) -> str:
    """Fetch a declared source beside its destination and move it in only once it matches.

    A partial transfer under the final name is the dangerous case: it looks like a
    finished download to every later stage. The temporary file is dot-prefixed so a
    crashed fetch does not appear as an input, and it is removed on any failure.
    """
    import urllib.request

    destination = Path(destination)
    temporary = destination.parent / ("." + destination.name + ".part")
    if temporary.exists():
        temporary.unlink()
    digest = hashlib.sha256()
    try:
        with urllib.request.urlopen(url, timeout=60) as response, temporary.open("wb") as handle:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                digest.update(chunk)
                handle.write(chunk)
    except Exception as exc:  # urllib raises a wide range of transport errors
        if temporary.exists():
            temporary.unlink()
        raise WorkflowError("download failed for %s: %s" % (url, exc))
    observed = digest.hexdigest()
    if observed != expected_sha256:
        temporary.unlink()
        raise WorkflowError("download of %s has sha256 %s, not the declared %s; the incomplete file "
                            "was removed" % (url, observed, expected_sha256))
    os.replace(str(temporary), str(destination))
    return observed


def outside_project(root: Path, value: Any) -> bool:
    """True when a declared input is written as a location outside the project.

    Deliberately lexical, and deliberately does not follow links: `Data/Matrix`
    pointing at a directory on shared storage is exactly the supported way to reach
    data stored elsewhere, so resolving it would report every correctly linked
    input as external. What this catches is the input that was never brought under
    `Data/` at all -- an absolute path elsewhere, or a relative path that climbs out.
    """
    if value in (None, ""):
        return False
    raw = Path(os.path.expandvars(os.path.expanduser(str(value))))
    root = Path(root)
    if raw.is_absolute():
        return root != raw and root not in raw.parents
    return ".." in raw.parts


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
            include = as_bool(raw["include"])
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


# --- root Report.md run log --------------------------------------------------
#
# The run log is a canonical re-render, never a surgical insert. The automatic
# path (inspect_run.py) and the rebuild path (update_report.py) both parse the
# existing region into sections, overlay the run summaries on disk, sort the
# union, and re-emit the whole region through build_report_text. Rebuilding is
# therefore identical to incrementally appending by construction rather than by
# keeping two code paths in agreement. It also removes a whole class of damage
# instead of repairing it: half markers, duplicated sections, and an END marker
# appearing before its START all converge to the same canonical text.

RUNLOG_PROLOGUE = (
    "<!-- 本区域由 tools/inspect_run.py 与 tools/update_report.py 生成并维护，请勿手工编辑；"
    "人类注记请写在本区域之上。 run-log format %d -->" % RUNLOG_FORMAT,
    "<!-- Generated and maintained by tools/inspect_run.py and tools/update_report.py; "
    "do not edit inside this region. Put human notes above it. -->",
)
RUNLOG_EMPTY = "_No run has been recorded yet. / 尚无运行记录。_"
RUNLOG_NOTE = (
    "_This file accumulates one concise record per analysis run. Re-inspecting a run "
    "updates its record in place; the section above is human text and is never rewritten. / "
    "本文件为每次分析运行累积一条精简记录；重跑同一 run 会就地更新其记录，区域上方为人类文本、不会被改写。_"
)
LEGACY_HEADING = "## Legacy run state / 旧运行状态（迁移保留，可手工删除）"
LEGACY_NOTE = (
    "_Recorded by an earlier release as a single overwritten block; kept verbatim. / "
    "由早前版本以单一覆写区块记录，原样保留。_"
)


def report_run_key(run_id: Any) -> str:
    """A marker key that can never break out of its HTML comment.

    plan_workflow only rejects "", "/", ".", and ".." for a run id, so ids such
    as `a-->b`, `a b`, or non-ASCII names reach the report. Such a key would
    terminate the comment and destroy marker parsing permanently, so anything
    outside the safe alphabet is replaced by a digest. The real run id is still
    shown verbatim in the record body.
    """
    text = str(run_id)
    if RUN_KEY_SAFE_RE.match(text):
        return text
    return "h" + hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _parse_checked_at(value: Any) -> Optional[datetime]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    # Mixed aware/naive stamps cannot be compared; pin naive stamps to UTC.
    return parsed if parsed.tzinfo is not None else parsed.replace(tzinfo=timezone.utc)


def run_record_lines(summary: Dict[str, Any], stage: Optional[str] = None) -> List[str]:
    """The concise per-run record: outcome, counts, and where the detail lives.

    With `stage`, only that stage's tasks are counted and the verdict is the
    stage's own. An empty list means the run never touched the stage, which is
    how a stage report leaves unrelated runs out instead of padding the log with
    records that say nothing about it.
    """
    run_id = str(summary.get("run_id") or "")
    status = str(summary.get("status") or "unknown")
    tasks = [row for row in (summary.get("tasks") or []) if isinstance(row, dict)]
    if stage is not None:
        # The summary written by inspect_run.py names the task `task`; older
        # summaries and in-memory rows use `id`. Both are read so a stage log
        # never silently drops a run it should have recorded.
        tasks = [row for row in tasks if stage_for_task(row.get("task") or row.get("id")) == stage]
        if not tasks:
            return []
        status = _stage_status(tasks)
    total = len(tasks)
    # Completion is proven by validated evidence, not by a scheduler state:
    # a not_submitted task and a marker-only task are both unfinished.
    complete = sum(1 for row in tasks if row.get("evidence_valid"))
    failed = sum(1 for row in tasks if str(row.get("state") or "") in TERMINAL_BAD_STATES)
    unfinished = max(0, total - complete - failed)
    signature = summary.get("input_signature")
    lines = [
        "### `%s` — %s" % (run_id, status),
        "",
        "- Status / 状态：`%s`" % status,
        "- Checked / 检查时间：`%s`" % str(summary.get("checked_at") or "unknown"),
        "- Tasks / 任务：%d/%d complete, %d failed, %d unfinished"
        % (complete, total, failed, unfinished),
        "- Input signature / 输入签名：`%s`" % (signature if signature else "unknown"),
        "- Evidence / 证据：`.workflow/runs/%s/run_summary.json`" % run_id,
    ]
    return lines


def render_run_section(key: str, lines: Iterable[str]) -> str:
    return "\n".join(
        ["<!-- SCMO-RUN:%s:START -->" % key]
        + list(lines)
        + ["<!-- SCMO-RUN:%s:END -->" % key]
    )


def _section(key: str, lines: List[str], index: int) -> Dict[str, Any]:
    body = list(lines)
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    run_id, when = None, None
    for line in body:
        if run_id is None:
            match = RUN_HEADING_RE.match(line)
            if match:
                run_id = match.group(1)
        if when is None:
            match = RUN_CHECKED_RE.match(line)
            if match:
                when = _parse_checked_at(match.group(1))
    return {"key": key, "lines": body, "run_id": run_id, "when": when, "index": index}


def _section_from_summary(summary: Dict[str, Any], index: int,
                          stage: Optional[str] = None) -> Optional[Dict[str, Any]]:
    lines = run_record_lines(summary, stage)
    if not lines:
        return None
    return {
        "key": report_run_key(summary.get("run_id") or ""), "lines": lines,
        "run_id": str(summary.get("run_id") or ""),
        "when": _parse_checked_at(summary.get("checked_at")), "index": index,
    }


def _section_sort_key(section: Dict[str, Any]) -> Tuple[int, float, int, str]:
    when = section.get("when")
    if when is None:
        # Hand-written or unparseable sections keep their file order, after the dated ones.
        return (1, 0.0, int(section.get("index") or 0), str(section.get("key") or ""))
    return (0, when.timestamp(), 0, str(section.get("key") or ""))


def parse_run_sections(region_text: str) -> Dict[str, Dict[str, Any]]:
    """Parse the per-run sections out of a run-log region.

    Deliberately total: every malformed arrangement yields a usable map instead
    of an exception, because the re-render repairs whatever it parsed.
    """
    sections: Dict[str, Dict[str, Any]] = {}
    open_key: Optional[str] = None
    buffer: List[str] = []
    index = 0
    for line in (region_text or "").splitlines():
        match = RUN_MARKER_RE.match(line)
        if not match:
            if open_key is not None:
                buffer.append(line)
            continue
        key, kind = match.group(1), match.group(2)
        if kind == "START":
            # An unterminated section is closed at the next START and re-emitted
            # with both markers, so the malformed state heals in one pass.
            if open_key is not None:
                sections[open_key] = _section(open_key, buffer, index)
                index += 1
            open_key, buffer = key, []
        elif open_key is not None and key == open_key:
            sections[open_key] = _section(open_key, buffer, index)
            index += 1
            open_key, buffer = None, []
        # A stray END with no matching open section is dropped.
    if open_key is not None:
        sections[open_key] = _section(open_key, buffer, index)
    return sections


def render_run_log_region(sections: Sequence[Dict[str, Any]]) -> str:
    lines = [RUNLOG_START] + list(RUNLOG_PROLOGUE)
    ordered = sorted(sections, key=_section_sort_key)
    if not ordered:
        lines += ["", RUNLOG_EMPTY]
    for section in ordered:
        lines += ["", render_run_section(str(section["key"]), section["lines"])]
    lines.append(RUNLOG_END)
    return "\n".join(lines)


def _migrate_legacy(text: str) -> str:
    """Replace the superseded single-block markers with a plain heading.

    Purely additive and idempotent: the legacy body stays exactly where it was,
    and once the two marker lines are gone this is a no-op on every later run.
    """
    out: List[str] = []
    migrated = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == LEGACY_RUN_STATE_START:
            out += [LEGACY_HEADING, "", LEGACY_NOTE, ""]
            migrated = True
            continue
        if stripped == LEGACY_RUN_STATE_END:
            migrated = True
            continue
        out.append(line)
    if not migrated:
        return text
    return "\n".join(out) + "\n"


def _drop_orphan_region_markers(text: str) -> str:
    """Remove region marker lines that have no partner.

    Reached only when no well-formed region exists, so every remaining region
    marker line is unpaired by definition and is tool-owned residue from an
    interrupted write. Dropping it here converges in one pass; leaving it would
    take a second refresh, because on the next read the orphan would be found
    first and swallow the freshly written region as body text.
    """
    kept = [line for line in text.splitlines() if not RUNLOG_MARKER_RE.match(line)]
    return "\n".join(kept)


def _split_region(text: str) -> Tuple[str, str, str]:
    """Split into (before the region, the region's body, after the region).

    The middle element is the region's inner text, with the markers themselves
    removed. Passing the markers on to parse_run_sections would let the closing
    RUNLOG_END be swept into the body of a section left open by a malformed
    write, and the re-render would then emit a second RUNLOG_END.
    """
    start = text.find(RUNLOG_START)
    if start == -1:
        return _drop_orphan_region_markers(text), "", ""
    # Searching for END only after START makes "END before START" impossible by
    # construction; testing the two markers independently is how the previous
    # implementation could slice a negative-length range and corrupt the file.
    end = text.find(RUNLOG_END, start + len(RUNLOG_START))
    if end == -1:
        return _drop_orphan_region_markers(text), "", ""
    return (
        text[:start],
        text[start + len(RUNLOG_START):end],
        text[end + len(RUNLOG_END):],
    )


def build_report_text(existing_text: str, summaries: Sequence[Dict[str, Any]],
                      stage: Optional[str] = None) -> str:
    """The single renderer for a Report's run log. Pure: no filesystem access.

    `stage` selects a stage view: the same summaries, filtered and counted per
    stage. The root Report and every stage report therefore render through one
    implementation, so a stage log is append-equal to the root log by
    construction rather than by two code paths agreeing.
    """
    migrated = _migrate_legacy(existing_text or "")
    # prefix and suffix are outside the region and belong to the human; only the
    # region's body is tool-owned and is rebuilt from scratch below.
    prefix, body, suffix = _split_region(migrated)
    # Sections already in the file survive a pruned or archived run directory;
    # the summaries on disk win on collision because they were just written. The
    # union is what makes the log accumulate rather than derive from disk.
    sections = parse_run_sections(body)
    for index, summary in enumerate(summaries or []):
        section = _section_from_summary(summary, index, stage)
        if section is not None:
            sections[section["key"]] = section
    rendered = render_run_log_region(list(sections.values()))
    # Blank lines are the one place a re-render could drift: the tail must lose
    # its trailing newlines too, or the next read of this file sees a different
    # junction and normalises to different bytes.
    head, tail = prefix.rstrip("\n"), suffix.strip("\n")
    parts = [part for part in (head, rendered, tail) if part]
    return "\n\n".join(parts) + "\n"


def collect_run_summaries(root: Path) -> Tuple[List[Dict[str, Any]], List[str]]:
    runs_dir = Path(root).resolve() / ".workflow" / "runs"
    summaries: List[Dict[str, Any]] = []
    warnings: List[str] = []
    if not runs_dir.is_dir():
        return summaries, warnings
    for path in sorted(runs_dir.glob("*/run_summary.json")):
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            warnings.append("%s is unreadable: %s" % (path, exc))
            continue
        if not isinstance(value, dict):
            warnings.append("%s is not a JSON object" % path)
            continue
        if not value.get("run_id"):
            value["run_id"] = path.parent.name
        summaries.append(value)
    return summaries, warnings


def _default_report_text(root: Path) -> str:
    project_id = root.name
    config_path = root / "config" / "project.yaml"
    if config_path.is_file():
        try:
            project_id = str((load_structured(config_path).get("project") or {}).get("id") or project_id)
        except (WorkflowError, OSError, ValueError):
            pass
    return (
        "# %s analysis report / 分析报告\n\n%s\n\n%s\n\n%s\n"
        % (project_id, RUNLOG_NOTE, RUNLOG_TITLE, render_run_log_region([]))
    )


def refresh_run_log(root: Path) -> Dict[str, Any]:
    """Re-render the root Report.md run-log region from the project's run summaries."""
    root = Path(root).resolve()
    report = root / "Report.md"
    existing = report.read_text(encoding="utf-8") if report.is_file() else _default_report_text(root)
    summaries, warnings = collect_run_summaries(root)
    text = build_report_text(existing, summaries)
    created = not report.is_file()
    if created or text != existing:
        write_text_atomic(report, text)
    return {
        "path": str(report),
        "records": len(re.findall(r"<!-- SCMO-RUN:[^\s<>]+:START -->", text)),
        "warnings": warnings,
        "written": created or text != existing,
    }


def refresh_stage_run_logs(root: Path) -> Dict[str, Any]:
    """Re-render every stage Report.md that the project already has.

    A stage report is a view of the same `.workflow/runs/*/run_summary.json` the
    root log reads. Reports for stages this project never generated a document
    for are skipped rather than created, because the human half of each file is
    that stage's contract and cannot be reconstructed from summaries.
    """
    root = Path(root).resolve()
    summaries, warnings = collect_run_summaries(root)
    stages: List[Dict[str, Any]] = []
    for stage, relative in STAGE_REPORTS:
        report = root / relative
        if not report.is_file():
            continue
        try:
            existing = report.read_text(encoding="utf-8")
        except OSError as exc:
            warnings.append("%s is unreadable: %s" % (report, exc))
            continue
        text = build_report_text(existing, summaries, stage=stage)
        if text != existing:
            write_text_atomic(report, text)
        stages.append({
            "stage": stage,
            "path": str(report),
            "records": len(re.findall(r"<!-- SCMO-RUN:[^\s<>]+:START -->", text)),
            "written": text != existing,
        })
    return {"stages": stages, "warnings": warnings}
