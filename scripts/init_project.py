#!/usr/bin/env python3
"""Create a standalone project from the bundled multiomics template."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence

from _common import (
    DATA_ROOT_NAME, RUNLOG_NOTE, RUNLOG_START, RUNLOG_TITLE, SAMPLE_COLUMNS, STAGE_CONTEXT_END,
    STAGE_CONTEXT_START, STAGE_REPORTS, WorkflowError, apply_data_sources, as_bool, data_sources,
    git_commit, load_structured, outside_project, plan_data_sources, render_run_log_region,
    resolve_path, write_json,
)


SKILL_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = SKILL_ROOT / "assets" / "project-template"
RUNTIME_TOOLS = (
    "_common.py", "bootstrap_environments.py", "validate_project.py", "plan_workflow.py", "inspect_resources.py",
    "submit_workflow.py", "inspect_run.py", "update_report.py", "record_annotation_review.py", "link_data.py",
)

# Reference files the skill ships so a methylation project needs no reference
# hunting before its first validation. They are copied into every generated
# project under Supplementary/ and are used only when the intake declares nothing
# for that key. The recorded digest is checked against the bundled file at
# generation time, because setting `blacklist_md5` from a file without looking at
# it would make the check that guards it unconditional instead of fail-closed.
BUNDLED_REFERENCES = (
    {"key": "chrom_sizes", "path": "Supplementary/hg38.canonical.chrom.sizes", "md5": None, "md5_key": None},
    {"key": "blacklist", "path": "Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz",
     "md5": "393688b4f06c9ce26165d47433dd8c37", "md5_key": "blacklist_md5"},
)
BUNDLED_GENOME_ID = "GRCh38"
HG38_GENOME_IDS = frozenset({"grch38", "hg38", "grch38.p14"})
HUMAN_ORGANISMS = frozenset({"human", "homo sapiens"})

DEFAULT_ANALYSIS = {
    # Every value mirrors the authoritative workflow notebook
    # (Scripts/Scanpy/Notebooks/scanpy_workflow.ipynb); the notebook documents each step.
    "scanpy": {
        "notebook": {
            "run_kind": "candidate", "iteration_id": None, "analysis_confirmed": False,
            "overwrite_data_outputs": True, "timeout_seconds": None, "output_stem": None,
        },
        "qc": {"min_genes": 200, "max_genes": 6000, "min_counts": 500,
               "max_mt_percent": 5.0, "min_cells_per_gene": 3},
        "scrublet": {"min_cells": 100, "doublet_rate_per_1000": 0.004},
        "target_sum": 10000,
        "hvg": {"n_top_genes": 2000, "batch_key": "cohort"},
        "pca": {"n_comps": 50, "svd_solver": "arpack"},
        "regress_covariates": False,
        "neighbors": {"n_pcs": 30, "n_neighbors": 15},
        "harmony": {"max_iter": 20, "sigma": 0.1},
        "umap": {"min_dist": 0.5, "spread": 1.0},
        "leiden": {"resolution": 0.8},
        "candidate": {"min_score_margin": 0.20, "top_markers_per_cluster": 5},
        "seed": 0,
        "sample_palette": ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"],
    },
    "methscan": {
        "min_sites": 300000, "min_meth_percent": 50, "max_meth_percent": 100,
        "smooth_bandwidth": 1000, "scan_bandwidth": 2000, "scan_stepsize": 100,
        "vmr_thresholds": [0.01, 0.02, 0.05], "min_cells": 6,
        "dmr": {"raw_p": 0.01, "min_abs_diff": 0.25, "top_per_cell_type": 200},
    },
    "allcools": {"bin_size": 5000, "mc_context": "CGN", "blacklist_fraction": 0.2},
    "methylvi": {
        "feature_targets": [10000, 30000], "epochs": 500, "batch_size": 32,
        "seed": 0, "validation_fraction": 0.1,
        "supervised_umap_weights": [0.2, 0.5, 0.7, 0.9],
    },
    "task_commands": {},
}

DEFAULT_PROFILES = {
    "default": {"floor": {"cpus": 1, "memory": "4G"}, "target": {"cpus": 4, "memory": "16G", "time": "1-00:00:00"}, "ceiling": {"cpus": 16, "memory": "64G"}},
    "serial": {"floor": {"cpus": 1, "memory": "8G"}, "target": {"cpus": 4, "memory": "16G", "time": "2-00:00:00"}, "ceiling": {"cpus": 8, "memory": "32G"}},
    "scanpy": {"floor": {"cpus": 4, "memory": "16G"}, "target": {"cpus": 16, "memory": "64G", "time": "2-00:00:00"}, "ceiling": {"cpus": 32, "memory": "128G"}},
    "io_builder": {"floor": {"cpus": 4, "memory": "16G"}, "target": {"cpus": 24, "memory": "64G", "time": "5-00:00:00"}, "ceiling": {"cpus": 56, "memory": "256G"}},
    "methscan_branch": {"floor": {"cpus": 4, "memory": "32G"}, "target": {"cpus": 18, "memory": "80G", "time": "2-00:00:00"}, "ceiling": {"cpus": 32, "memory": "160G"}},
    # The `dmr` floor must not fall below the parallelism the meth-diff stage declares for itself:
    # 07_methdiff_celltype.py fails the task when SLURM_CPUS_PER_TASK is less than jobs*threads
    # (2*4 by default), so a floor of 4 let a busy node hand back fewer CPUs than the stage needs.
    # Raising SCMO_METHDIFF_JOBS or SCMO_METHDIFF_THREADS means raising this floor to their product.
    "dmr": {"floor": {"cpus": 8, "memory": "16G"}, "target": {"cpus": 8, "memory": "32G", "time": "5-00:00:00"}, "ceiling": {"cpus": 32, "memory": "96G"}},
    "dmr_prepare": {"floor": {"cpus": 2, "memory": "16G"}, "target": {"cpus": 8, "memory": "32G", "time": "1-00:00:00"}, "ceiling": {"cpus": 16, "memory": "64G"}},
    "feature_builder": {"floor": {"cpus": 8, "memory": "32G"}, "target": {"cpus": 32, "memory": "90G", "time": "7-00:00:00"}, "ceiling": {"cpus": 64, "memory": "256G"}},
    "trainer": {"floor": {"cpus": 8, "memory": "32G"}, "target": {"cpus": 32, "memory": "90G", "time": "7-00:00:00"}, "ceiling": {"cpus": 64, "memory": "256G"}},
    "plot": {"floor": {"cpus": 2, "memory": "8G"}, "target": {"cpus": 8, "memory": "32G", "time": "1-00:00:00"}, "ceiling": {"cpus": 16, "memory": "64G"}},
    "summary": {"floor": {"cpus": 1, "memory": "2G"}, "target": {"cpus": 2, "memory": "8G", "time": "02:00:00"}, "ceiling": {"cpus": 4, "memory": "16G"}},
}


def merge(base, override):
    result = json.loads(json.dumps(base))
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def bundled_references(intake_references, organism, samples):
    """Fill in the reference files the skill ships and report which ones were used.

    Three conditions have to hold before inheriting a bundled human reference is
    honest, and each of them is a way this could go quietly wrong rather than
    loudly: a project with no ALLC input never reads them (dead weight that still
    lands in every evidence record), a non-human project would filter the wrong
    regions while looking perfectly healthy, and a GRCh37 study would inherit
    GRCh38 coordinates. A project that fails one of them is told why and keeps the
    preflight failure that names the reference it must declare.
    """
    references = dict(intake_references or {})
    notes: List[Dict[str, Any]] = []
    warnings: List[str] = []
    methylation = any(
        isinstance(sample, dict) and as_bool(sample.get("include")) and str(sample.get("allc_root") or "").strip()
        for sample in samples or []
    )
    if not methylation:
        return references, notes, warnings
    # Checked in lowercase, reported as written: the warning quotes what the intake
    # actually says, so it cannot be mistaken for a value the reader never wrote.
    organism = str(organism or "").strip()
    genome = str(references.get("genome_id") or "").strip()
    if organism and organism.lower() not in HUMAN_ORGANISMS:
        warnings.append(
            "no bundled genome references were applied: this project declares organism %r and the bundled "
            "references are human hg38. Declare references.chrom_sizes and references.blacklist in the intake."
            % organism)
        return references, notes, warnings
    if genome and genome.lower() not in HG38_GENOME_IDS:
        warnings.append(
            "no bundled genome references were applied: this project declares genome_id %r and the bundled "
            "references are hg38. Declare references.chrom_sizes and references.blacklist in the intake."
            % genome)
        return references, notes, warnings
    for entry in BUNDLED_REFERENCES:
        if references.get(entry["key"]):
            continue
        bundled = TEMPLATE / entry["path"]
        if not bundled.is_file():
            raise WorkflowError("the skill is missing the bundled reference %s" % bundled)
        if entry["md5"] and hashlib.md5(bundled.read_bytes()).hexdigest() != entry["md5"]:
            raise WorkflowError(
                "bundled reference %s does not match its recorded digest %s; the skill install is corrupt"
                % (bundled, entry["md5"]))
        references[entry["key"]] = entry["path"]
        if entry["md5"] and entry["md5_key"]:
            references[entry["md5_key"]] = entry["md5"]
        notes.append(entry)
    if notes and not references.get("genome_id"):
        # Written so the config is self-consistent with the hg38 files it now
        # points at, and recorded as a note rather than silently: the README
        # renders this key as a project fact, so the Report has to say where it
        # came from. Nothing reads genome_id today; a reader does.
        references["genome_id"] = BUNDLED_GENOME_ID
        notes.append({"inferred_genome": BUNDLED_GENOME_ID})
    return references, notes, warnings


# --- root README -------------------------------------------------------------
#
# The README is written once, from the intake and the resolved configuration, and
# is never rewritten with runtime state: its job is to be the stable contract a
# reader can trust months later. Everything that changes as analysis proceeds
# belongs in Report.md. Both halves of the document are rendered from the same
# row tuples so the two languages cannot drift apart.
README_PARAM_ROWS = (
    ("质控：最少基因数", "QC: min genes", "scanpy.qc.min_genes"),
    ("质控：最多基因数", "QC: max genes", "scanpy.qc.max_genes"),
    ("质控：最少 counts", "QC: min counts", "scanpy.qc.min_counts"),
    ("质控：线粒体占比上限", "QC: max mitochondrial fraction", "scanpy.qc.max_mt_percent"),
    ("高变基因数", "Highly variable genes", "scanpy.hvg.n_top_genes"),
    ("高变基因批次键", "HVG batch key", "scanpy.hvg.batch_key"),
    ("PCA 维数", "PCA components", "scanpy.pca.n_comps"),
    ("近邻：PC 数", "Neighbors: PCs", "scanpy.neighbors.n_pcs"),
    ("近邻：邻居数", "Neighbors: k", "scanpy.neighbors.n_neighbors"),
    ("Harmony 迭代数 / sigma", "Harmony iterations / sigma", "scanpy.harmony.max_iter"),
    ("Leiden 分辨率", "Leiden resolution", "scanpy.leiden.resolution"),
    ("目标 counts 归一化", "Target sum", "scanpy.target_sum"),
    ("随机种子", "Seed", "scanpy.seed"),
    ("MethSCAn：最少位点", "MethSCAn: min sites", "methscan.min_sites"),
    ("MethSCAn：VMR 阈值", "MethSCAn: VMR thresholds", "methscan.vmr_thresholds"),
    ("MethSCAn：最少细胞", "MethSCAn: min cells", "methscan.min_cells"),
    ("MethSCAn DMR：p 值 / 最小差异", "MethSCAn DMR: p / min diff", "methscan.dmr.raw_p"),
    ("ALLCools bin 大小", "ALLCools bin size", "allcools.bin_size"),
    ("ALLCools 甲基化上下文", "ALLCools context", "allcools.mc_context"),
    ("MethylVI 特征目标", "MethylVI feature targets", "methylvi.feature_targets"),
    ("MethylVI epochs / batch size", "MethylVI epochs / batch size", "methylvi.epochs"),
    ("MethylVI 验证集比例", "MethylVI validation fraction", "methylvi.validation_fraction"),
)


def _dotted(config, path: str):
    """Read a dotted path, or None when any level is absent or not a mapping."""
    node = config
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def _cell(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


def _profile_cell(profile, level: str) -> str:
    tier = (profile or {}).get(level) or {}
    parts = [str(tier.get("cpus", "?")) + "cpu", str(tier.get("memory", "?"))]
    if level == "target" and tier.get("time"):
        parts.append(str(tier["time"]))
    return " / ".join(parts)


def _sample_rows(samples) -> List[List[str]]:
    rows: List[List[str]] = []
    for sample in samples or []:
        if not isinstance(sample, dict):
            continue
        rna = str(sample.get("rna_path") or "").strip()
        allc = str(sample.get("allc_root") or "").strip()
        modality = "RNA+ALLC" if rna and allc else "RNA" if rna else "ALLC" if allc else "—"
        rows.append([
            _cell(sample.get("sample_id")),
            _cell(sample.get("condition")),
            _cell(sample.get("batch")),
            "yes" if as_bool(sample.get("include")) else "no",
            modality,
        ])
    return rows


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return lines


def _readme_half(analysis_cfg, scheduler_cfg, samples, environments, paths, lang, identity, data_rows, data_outside) -> List[str]:
    """One language half, rendered from the shared bilingual row tables.

    `lang` selects which side of every (zh, en) pair to show, so the halves stay
    in step without duplicating any of the data or the section structure.
    """
    scheduler = scheduler_cfg or {}
    rows = _sample_rows(samples)
    params = []
    for zh, en, key in README_PARAM_ROWS:
        value = _dotted(analysis_cfg or {}, key)
        if value is not None:
            params.append([zh if lang == "zh" else en, "`%s`" % key, "`%s`" % _cell(value)])
    profiles = scheduler.get("profiles") or {}
    profile_rows = [
        [name, _profile_cell(profiles[name], "floor"), _profile_cell(profiles[name], "target"),
         _profile_cell(profiles[name], "ceiling")]
        for name in sorted(profiles) if isinstance(profiles.get(name), dict)
    ]
    lines = list(paths["intro"])
    lines += ["", paths["identity_heading"], ""] + _table(paths["identity_header"], identity)
    lines += ["", paths["samples_heading"], ""]
    if rows:
        lines += _table(paths["sample_header"], rows)
        lines += ["", paths["samples_pointer"]]
    else:
        lines += [paths["no_samples"]]
    lines += ["", paths["data_heading"], ""]
    if data_rows:
        lines += _table(paths["data_header"],
                        [[row["name"], "`%s`" % row["origin"], paths["data_states"][row["state"]]]
                         for row in data_rows])
    else:
        lines += [paths["no_data"]]
    lines += ["", paths["data_pointer"]]
    if data_outside:
        lines += ["", paths["data_outside"] % ", ".join("`%s`" % item for item in data_outside)]
    lines += ["", paths["params_heading"], ""]
    lines += _table(paths["param_header"], params) if params else [paths["no_params"]]
    lines += ["", paths["params_pointer"], "", paths["env_heading"], ""]
    if environments:
        env_rows = [
            [_cell(row.get("stage")), _cell(row.get("python")), _cell(row.get("required"))]
            for row in environments if isinstance(row, dict)
        ]
        lines += _table(paths["env_header"], env_rows)
    else:
        lines += [paths["no_env"]]
    lines += ["", paths["sched_heading"], ""]
    lines += _table(paths["sched_header"], [
        [paths["backend_label"], "`%s`" % _cell(scheduler.get("backend"))],
        [paths["account_label"], "`%s`" % _cell(scheduler.get("account")) if scheduler.get("account") else paths["account_unset"]],
        [paths["partitions_label"], "`%s`" % _cell(scheduler.get("partitions")) if scheduler.get("partitions") else paths["partitions_unset"]],
        [paths["parallel_label"], "`%s`" % _cell(scheduler.get("max_parallel"))],
    ])
    lines += ["", paths["profiles_pointer"], ""]
    lines += _table(paths["profile_header"], profile_rows) if profile_rows else [paths["no_profiles"]]
    lines += ["", paths["run_heading"], ""]
    lines += ["- " + item for item in paths["run_items"]]
    lines += ["", paths["frozen"]]
    return lines


def bilingual_readme(project_cfg, analysis_cfg, scheduler_cfg, samples, environments,
                     data_rows=None, data_outside=None) -> str:
    """Render the root README from the project's own facts.

    Deliberately `%`-free and `.format`-free assembly: the document quotes JSON
    literals, and `str.format` would require doubling every brace in them.
    Every lookup tolerates a missing key, because this runs for schema-v1
    intakes and sparse fixtures too.
    """
    project_id = str((project_cfg.get("project") or {}).get("id") or "project")
    project = project_cfg.get("project") or {}
    references = project_cfg.get("references") or {}
    zh_paths = {
        "intro": [
            "本项目由 `single-cell-multiomics-analysis` 生成。所有路径、样本、参考基因组、环境和调度资源均由 `config/` 显式定义。缺失环境先通过 `tools/bootstrap_environments.py --execute` 隔离创建并验证，再执行 quick/full 校验；正式任务只能通过 `tools/submit_workflow.py` 提交，结果写入 `Results/runs/<run-id>/`。专家覆写任务必须生成 `<task_dir>/task_outputs.json`，其中 `artifacts` 是非空且实际存在的绝对路径列表。",
        ],
        "identity_heading": "### 项目",
        "identity_header": ["项", "值"],
        "samples_heading": "### 样本",
        "sample_header": ["样本 ID", "条件", "批次", "纳入分析", "模态"],
        "samples_pointer": "`config/samples.tsv` 是完整清单（含 RNA 路径、ALLC 根目录、细胞 ID 规则等全部列）；上表只列出识别项目所需的字段。",
        "no_samples": "intake 未声明样本；请先补充 `config/samples.tsv`。",
        "data_heading": "### 输入数据",
        "data_header": ["`Data/` 条目", "来源", "状态"],
        "data_states": {"in_place": "已就位", "dangling": "已链接；生成主机上不可见",
                        "pending": "待下载", "unlinked": "已声明，未链接（`--no-link`）"},
        "data_pointer": "分析输入一律从项目内的 `Data/` 读取：在 intake 的 `data.sources` 中声明来源，生成时软链接为 `Data/<name>`；`url` 来源改由 `tools/link_data.py --execute` 下载并校验 sha256。`config/samples.tsv` 中的路径据此写成 `Data/...`。",
        "no_data": "intake 未声明 `data.sources`；`config/samples.tsv` 中的路径即各输入的实际位置。",
        "data_outside": "以下样本读取项目外的路径而非经由 `Data/`：%s。分析照常读取，但项目因此不再自包含；建议在 intake 的 `data.sources` 中声明来源，使其落到 `Data/` 下。",
        "params_heading": "### 关键分析参数",
        "param_header": ["参数", "配置键", "值"],
        "params_pointer": "上表是解读结果所需的关键参数，不是全部：完整且权威的参数在 `config/analysis.yaml`。",
        "no_params": "`config/analysis.yaml` 中无可显示参数。",
        "env_heading": "### 隔离环境",
        "env_header": ["阶段", "Python", "必需"],
        "no_env": "生成时未声明隔离环境。请运行 `tools/bootstrap_environments.py --execute` 创建并验证后再做正式提交。",
        "sched_heading": "### 调度资源",
        "sched_header": ["项", "值"],
        "backend_label": "后端", "account_label": "账户",
        "partitions_label": "分区", "parallel_label": "最大并行",
        "account_unset": "未设置（按站点填写 `config/scheduler.yaml`）",
        "partitions_unset": "未设置（按站点填写 `config/scheduler.yaml`）",
        "profiles_pointer": "命名档位定义在 `config/scheduler.yaml`；档位是唯一资源来源，阶段脚本本身不含站点相关的 `#SBATCH` 指令。",
        "profile_header": ["档位", "下限", "目标", "上限"],
        "no_profiles": "未定义命名档位。",
        "run_heading": "### 运行",
        "run_items": [
            "环境：`python tools/bootstrap_environments.py --execute`",
            "校验：quick 校验用于快速体检，full 校验是正式提交的前置；结论记录在 `.workflow/validations/`。",
            "计划：`python tools/plan_workflow.py --project . --routes auto`",
            "提交：`python tools/submit_workflow.py --project . --run-id <run-id>`",
            "检查：`python tools/inspect_run.py --project . --run-id <run-id>`，同时把该 run 的精简记录写入 `Report.md`。",
        ],
        "frozen": "**本文件由 `init_project.py` 依据 intake 与配置写入一次，此后永不写入运行时状态**（不反映进度、不记录结果）。运行台账在 `Report.md` 中按 run 累积；模块说明见 `Scripts/<Module>/README.md`。",
    }
    en_paths = {
        "intro": [
            "This project was generated by `single-cell-multiomics-analysis`. Paths, samples, genome resources, environments, and scheduler resources are declared under `config/`. Bootstrap missing isolated environments, run quick/full validation, and submit production jobs only through `tools/submit_workflow.py`. Results are isolated under `Results/runs/<run-id>/`. Expert overrides must create `<task_dir>/task_outputs.json` with a non-empty `artifacts` list of existing absolute paths.",
        ],
        "identity_heading": "### Project",
        "identity_header": ["Item", "Value"],
        "samples_heading": "### Samples",
        "sample_header": ["Sample ID", "Condition", "Batch", "Included", "Modality"],
        "samples_pointer": "`config/samples.tsv` is the complete manifest (RNA paths, ALLC roots, cell-id rules, and the remaining columns); the table above lists only the identifying fields.",
        "no_samples": "The intake declared no samples; add them to `config/samples.tsv`.",
        "data_heading": "### Input data",
        "data_header": ["`Data/` entry", "Origin", "State"],
        "data_states": {"in_place": "in place", "dangling": "linked; not visible from the generating host",
                        "pending": "pending download", "unlinked": "declared, not linked (`--no-link`)"},
        "data_pointer": "Analysis inputs are read only from the project's own `Data/`: declare their origin under `data.sources` in the intake and generation links it to `Data/<name>`, or fetch a `url` source with `tools/link_data.py --execute`, which verifies the sha256. Paths in `config/samples.tsv` are written as `Data/...` accordingly.",
        "no_data": "The intake declared no `data.sources`; the paths in `config/samples.tsv` are the inputs' actual locations.",
        "data_outside": "These samples read paths outside the project instead of through `Data/`: %s. The analysis still reads them, but the project is no longer self-contained; declare their origin in the intake's `data.sources` so they land under `Data/`.",
        "params_heading": "### Key analysis parameters",
        "param_header": ["Parameter", "Config key", "Value"],
        "params_pointer": "The table above holds the parameters needed to interpret results, not all of them: the complete, authoritative set is in `config/analysis.yaml`.",
        "no_params": "`config/analysis.yaml` exposes no displayable parameters.",
        "env_heading": "### Isolated environments",
        "env_header": ["Stage", "Python", "Required"],
        "no_env": "The intake declared no isolated environments. Run `tools/bootstrap_environments.py --execute` to create and verify them before submitting production work.",
        "sched_heading": "### Scheduler resources",
        "sched_header": ["Item", "Value"],
        "backend_label": "Backend", "account_label": "Account",
        "partitions_label": "Partitions", "parallel_label": "Max parallel",
        "account_unset": "unset (fill in `config/scheduler.yaml` for your site)",
        "partitions_unset": "unset (fill in `config/scheduler.yaml` for your site)",
        "profiles_pointer": "Named profiles are defined in `config/scheduler.yaml`; profiles are the only source of resources and stage scripts carry no site-specific `#SBATCH` directives.",
        "profile_header": ["Profile", "Floor", "Target", "Ceiling"],
        "no_profiles": "No named profiles are defined.",
        "run_heading": "### Running",
        "run_items": [
            "Environments: `python tools/bootstrap_environments.py --execute`",
            "Validation: quick validates cheaply, full validation gates a production submit; verdicts are recorded under `.workflow/validations/`.",
            "Plan: `python tools/plan_workflow.py --project . --routes auto`",
            "Submit: `python tools/submit_workflow.py --project . --run-id <run-id>`",
            "Inspect: `python tools/inspect_run.py --project . --run-id <run-id>`, which also writes this run's concise record into `Report.md`.",
        ],
        "frozen": "**This file is written once by `init_project.py` from the intake and configuration and is never rewritten with runtime state** (it reports no progress and no results). The run log accumulates per run in `Report.md`; module documentation is in `Scripts/<Module>/README.md`.",
    }
    # Identity rows are language-label pairs so both halves stay in step.
    identity_rows = [
        ("项目 ID", "Project ID", project_id),
        ("生物", "Organism", project.get("organism")),
        ("参考基因组", "Genome", project.get("genome_id")),
        ("项目根", "Project root", project.get("root")),
    ]
    for key in ("annotation", "chrom_sizes", "blacklist", "genome_fasta"):
        if references.get(key):
            identity_rows.append(("参考：" + key, "Reference: " + key, references.get(key)))
    lines = ["# " + project_id]
    for label, lang, paths in (("## 中文", "zh", zh_paths), ("## English", "en", en_paths)):
        identity = [
            [zh if lang == "zh" else en, "`%s`" % _cell(value)]
            for zh, en, value in identity_rows
        ]
        lines += ["", label]
        lines += [""] + _readme_half(
            analysis_cfg, scheduler_cfg, samples, environments, paths, lang, identity,
            data_rows or [], data_outside or [],
        )
    return "\n".join(lines) + "\n"


def bilingual_report(project_id: str, provenance: Dict[str, Any], data_rows=None,
                     data_outside=None, reference_notes=None, reference_warnings=None) -> str:
    """Render the root Report.md: the human plan, then the tool-owned run region.

    The region is emitted through the same renderer the tools use, so the file is
    canonical from its first byte and an inspection that changes nothing rewrites
    nothing. Everything above the region is written once, here, and records what
    the project was generated from -- above all which inputs and which bundled
    reference files it inherited, which nothing else in the project states.
    """
    planned = [
        "- Status / 状态: initialized",
        "- Generated / 生成时间: `%s`" % _cell((provenance or {}).get("generated_at")),
        "- Skill / 技能: `single-cell-multiomics-analysis` %s" % _cell((provenance or {}).get("skill_version")),
        "- Skill commit / 技能提交: `%s`" % _cell((provenance or {}).get("skill_git_commit")),
        "- Runtime evidence / 运行证据: `.workflow/runs/<run-id>/run_summary.json`",
        "- 本段位于运行记录区域之上，属人类文本，工具不会改写。",
        "- This section sits above the run-log region and is human text; the tools never rewrite it.",
    ]
    inputs = ["## Inputs / 输入数据", ""]
    if data_rows:
        inputs += ["| `Data/` entry | Origin / 来源 | State / 状态 |", "|---|---|---|"]
        inputs += ["| `%s` | `%s` | %s |" % (row["name"], row["origin"], row["state"]) for row in data_rows]
    else:
        inputs += ["- 未声明 `data.sources`；输入位置直接写在 `config/samples.tsv`。",
                   "- No `data.sources` declared; input locations are written directly in `config/samples.tsv`."]
    if data_outside:
        inputs += ["",
                   "- 以下样本读取项目外路径而非经由 `Data/`：%s。" % ", ".join("`%s`" % item for item in data_outside),
                   "- These samples read paths outside the project rather than through `Data/`: %s."
                   % ", ".join("`%s`" % item for item in data_outside)]
    references = ["## Bundled references / 内置参考默认值", ""]
    if reference_notes:
        references += [
            "- 以下参考文件取自技能内置资源（未在 intake 中声明）：",
            "- These reference files come from the skill's bundled resources (the intake declared none):",
        ]
        for entry in reference_notes:
            if entry.get("inferred_genome"):
                references += [
                    "  - 未声明基因组，而内置默认值基于 hg38，故 `references.genome_id` 记为 `%s`；"
                    "数据若不是 hg38，请连同参考文件一并覆盖。"
                    % entry["inferred_genome"],
                    "  - The intake declared no genome and the bundled defaults are hg38, so "
                    "`references.genome_id` records `%s`; override it along with the reference files if "
                    "the data is not hg38." % entry["inferred_genome"],
                ]
                continue
            digest = " (md5 `%s`)" % entry["md5"] if entry["md5"] else ""
            references += ["  - `references.%s` = `%s`%s" % (entry["key"], entry["path"], digest)]
        references += ["",
                       "- 这是生成时写入的事实记录；在 `config/project.yaml` 的 `references` 中声明自己的文件即可覆盖，"
                       "但那会改变输入签名，已记录的全量校验需要重做。",
                       "- This records what generation did. Declaring your own files under `references` in "
                       "`config/project.yaml` overrides it, at the cost of a new input signature and a repeat of "
                       "the recorded full validation."]
    if reference_warnings:
        references += [""] + ["- " + warning for warning in reference_warnings]
    if not reference_notes and not reference_warnings:
        references += ["- 未使用内置参考文件。", "- No bundled reference files are in use."]
    parts = [
        "# %s analysis report / 分析报告" % project_id,
        "## Planned / 计划\n\n" + "\n".join(planned),
        "\n".join(inputs),
        "\n".join(references),
        RUNLOG_NOTE,
        RUNLOG_TITLE,
        render_run_log_region([]),
    ]
    return "\n\n".join(parts) + "\n"


# --- packaged stage documents -------------------------------------------------
#
# The stage READMEs and Reports are bilingual contracts ported verbatim from the
# reference project. Their prose keeps that project's example sample names on
# purpose: the figures printed beside those names are that project's measured
# values, and rewriting them into statements about different data would fabricate
# evidence. What generation can honestly do is say which project this copy
# belongs to, so a reader never has to guess whether `CYL`/`ZCP` are examples or
# their own samples. This block is that statement, written once and never
# rewritten by a run, and the stage Reports additionally start with an empty run
# log that every later inspection fills in.

STAGE_RUNLOG_NOTE = (
    "_This file accumulates one concise record per run that touched this stage. Re-inspecting "
    "a run updates its record in place; everything above the region is human text and is never "
    "rewritten. The root `Report.md` holds the same records for every stage together. / "
    "本文件为该阶段每次运行累积一条精简记录；重跑同一 run 会就地更新其记录，区域上方为人类文本、"
    "不会被改写。根目录 `Report.md` 汇总全部阶段的相同记录。_"
)


def packaged_notebook(output: Path) -> str:
    """The stem of the notebook this project runs, from the packaged manifest.

    The manifest is the one place that declares which notebook is the entry, so
    reading it keeps the documents and the runner from drifting apart if the
    template notebook is ever renamed.
    """
    manifest = output / "Scripts" / "Scanpy" / "Notebooks" / "scanpy_template_manifest.json"
    try:
        configured = str(json.loads(manifest.read_text(encoding="utf-8")).get("template") or "")
    except (OSError, ValueError, TypeError):
        return "scanpy_workflow"
    return configured[:-6] if configured.endswith(".ipynb") else (configured or "scanpy_workflow")


def stage_context_block(project_id: str, organism: Any, genome_id: Any, labels: Sequence[str],
                        notebook: str) -> str:
    """The generated header naming what this project actually is."""
    facts = [
        "- 项目 / Project: `%s`" % project_id,
        "- 物种 / Organism: %s" % (_cell(organism) if organism else "—"),
        "- 基因组 / Genome: `%s`" % (genome_id if genome_id else "—"),
    ]
    if labels:
        facts.append("- 样本 / Samples (%d): %s"
                     % (len(labels), ", ".join("`%s`" % label for label in labels)))
    if notebook:
        facts.append("- Scanpy 入口 / Scanpy entry: `Notebooks/%s.ipynb`" % notebook)
    lines = [STAGE_CONTEXT_START, "> **本项目上下文 / Project context**", ">"]
    lines += ["> " + fact for fact in facts]
    lines.append(">")
    lines.append(
        "> 本段由 `init_project.py` 生成一次，不随运行改写。"
        + ("下文出现的 `CYL`/`ZCP` 等样本名属于参考项目，本项目的实际样本以上表为准。"
           if labels else "")
    )
    lines.append(
        "> Generated once by `init_project.py` and never rewritten by a run. "
        + ("Sample names such as `CYL`/`ZCP` below belong to the reference project; this "
           "project's own samples are the ones listed above."
           if labels else "Numbers and figures below describe the reference project; this "
           "project's own values are recorded in its run log.")
    )
    lines.append(STAGE_CONTEXT_END)
    return "\n".join(lines)


def adapt_stage_docs(output: Path, project_id: str, organism: Any, genome_id: Any,
                     labels: Sequence[str], notebook: str) -> List[str]:
    """Bind every packaged stage document to this project, once, at generation.

    Two edits, both confined to what generation actually knows: the context block
    above, and the `<notebook>` placeholder resolved to the notebook this project
    ships. Idempotent, so re-running generation over an existing tree is safe.

    Deliberately narrow. The stage documents are contracts ported from the
    reference project, and their worked numbers are that project's measurements;
    rewriting those into statements about this project would invent evidence. The
    banner says whose numbers they are, and the run log below says what actually
    happened here.
    """
    banner = stage_context_block(project_id, organism, genome_id, labels, notebook)
    # Only reports a run actually writes get a run-log region. A region elsewhere
    # would claim records nothing ever fills.
    logged = {relative for _, relative in STAGE_REPORTS}
    adapted: List[str] = []
    for path in sorted((output / "Scripts").rglob("*.md")):
        if path.name not in {"README.md", "Report.md"}:
            continue
        relative = path.relative_to(output).as_posix()
        text = path.read_text(encoding="utf-8")
        if STAGE_CONTEXT_START in text:
            continue
        lines = text.splitlines()
        # After the H1 when there is one, so the block reads as part of the header.
        index = 1 if lines and lines[0].startswith("# ") else 0
        text = "\n".join(lines[:index] + ["", banner] + lines[index:]).rstrip("\n") + "\n"
        if notebook:
            text = text.replace("Notebooks/<notebook>.ipynb", "Notebooks/%s.ipynb" % notebook)
        if relative in logged and RUNLOG_START not in text:
            # Seed an empty region so the first inspection appends into a placed
            # region instead of bolting one onto the end of the document.
            text = "\n".join([text.rstrip("\n"), "", STAGE_RUNLOG_NOTE, "", RUNLOG_TITLE, "",
                              render_run_log_region([])]) + "\n"
        path.write_text(text, encoding="utf-8")
        adapted.append(relative)
    return adapted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intake", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-link", action="store_true",
                        help="declare data.sources without creating the Data/ links they describe")
    args = parser.parse_args()
    intake = load_structured(args.intake)
    if intake.get("schema_version") not in {1, 2}:
        raise WorkflowError("intake schema_version must be 1 or 2")
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise WorkflowError("refusing non-empty output directory: %s" % output)
    output.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        str(TEMPLATE), str(output), dirs_exist_ok=True,
        ignore=shutil.ignore_patterns(
            "__pycache__", "*.pyc", "*.pyo", ".pytest_cache",
            ".ipynb_checkpoints", ".DS_Store", "Config",
        ),
    )
    tools_dir = output / "tools"
    tools_dir.mkdir(exist_ok=True)
    for name in RUNTIME_TOOLS:
        shutil.copy2(str(SKILL_ROOT / "scripts" / name), str(tools_dir / name))

    project_id = str(intake.get("project_id", "")).strip()
    if not project_id:
        raise WorkflowError("project_id is required")
    # A relative data.sources path is anchored to the intake file, so the same
    # intake links the same targets regardless of the directory it is run from.
    entries = data_sources(intake.get("data"), args.intake.resolve().parent)
    references, reference_notes, reference_warnings = bundled_references(
        intake.get("references"), intake.get("organism"), intake.get("samples"))
    project_cfg = {
        "schema_version": 2,
        "project": {
            "id": project_id, "root": str(output), "language": "zh-en",
            "organism": intake.get("organism", "unspecified"),
            "genome_id": references.get("genome_id", "unspecified"),
            "mitochondrial_prefixes": intake.get("mitochondrial_prefixes", ["MT-", "mt-"]),
            "ribosomal_prefixes": intake.get("ribosomal_prefixes", ["RPL", "RPS"]),
        },
        "paths": {"results": "Results", "logs": "Scripts/logs", "state": ".workflow"},
        "data": {"root": DATA_ROOT_NAME, "sources": entries},
        "references": references,
        "annotation": intake.get("annotation") or {
            "table": None, "profile": None, "review_status": "unreviewed",
            "cell_id_column": "cell_id", "cell_type_column": "cell_type",
        },
    }
    annotation = project_cfg["annotation"]
    if "path" in annotation and "table" not in annotation:
        annotation["table"] = annotation["path"]
    annotation.setdefault("profile", None)
    annotation.setdefault("review_status", "unreviewed")
    analysis_cfg = {"schema_version": 2, "analysis": merge(DEFAULT_ANALYSIS, intake.get("analysis") or {})}
    scheduler_cfg = merge({
        "schema_version": 2, "backend": "local", "account": None,
        "partitions": [], "allow_nodes": [], "exclude_nodes": [],
        "memory_headroom_mb": 4096, "max_parallel": 2, "validation_workers": 4,
        "limited_profiles": ["methscan_branch", "dmr", "feature_builder", "trainer"],
        "local": {"max_threads": 16, "max_memory": "64G"}, "profiles": DEFAULT_PROFILES,
    }, intake.get("scheduler") or {})
    scheduler_cfg["schema_version"] = 2
    write_json(output / "config" / "project.yaml", project_cfg)
    write_json(output / "config" / "analysis.yaml", analysis_cfg)
    write_json(output / "config" / "scheduler.yaml", scheduler_cfg)

    samples = intake.get("samples") or []
    with (output / "config" / "samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in samples:
            writer.writerow({key: sample.get(key, "") for key in SAMPLE_COLUMNS})
    env_rows = intake.get("environments") or []
    env_columns = ("stage", "python", "executable", "version_command", "required")
    with (output / "config" / "environments.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=env_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in env_rows:
            writer.writerow({key: row.get(key, "") for key in env_columns})
    # One provenance record feeds both the report header and the template lock,
    # so the two can never claim different skill versions.
    provenance = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "skill_version": (SKILL_ROOT / "VERSION").read_text().strip(),
        "skill_git_commit": git_commit(SKILL_ROOT),
    }
    # Links are created here because the intake's paths are visible on the host
    # that has just read them. Downloads are not: a url source can be arbitrarily
    # large, so it is recorded as pending and fetched by tools/link_data.py, which
    # says what it is about to transfer before it transfers it.
    for directory in (output / "Results" / "runs", output / "Supplementary",
                      output / "Scripts" / "logs", output / DATA_ROOT_NAME):
        directory.mkdir(parents=True, exist_ok=True)
    data_rows, data_outside = [], []
    if entries and not args.no_link:
        for result in apply_data_sources(output, plan_data_sources(output, entries)):
            name = result["name"]
            data_rows.append({
                "name": name, "origin": result["url"] or result["target"],
                "state": "pending" if result["outcome"] == "pending" else
                         "dangling" if result["outcome"] in {"linked-unresolved", "unresolved"} else "in_place",
            })
            if result["outcome"] == "pending":
                print("NOTE: %s is declared as a download and is not fetched yet; run "
                      "`python tools/link_data.py --project . --execute`" % name, file=sys.stderr)
            elif result["outcome"] in {"linked-unresolved", "unresolved"}:
                print("WARNING: Data/%s points at %s, which does not exist from this host; if that path is "
                      "only absent here, validate from the execution context" % (name, result["target"]),
                      file=sys.stderr)
    elif entries:
        # --no-link: declared but deliberately not created here, which is a third
        # state. Calling it "pending" would read as a download that has not run.
        data_rows = [{"name": entry["name"], "origin": entry["url"] or entry["target"], "state": "unlinked"}
                     for entry in entries]
    for sample in samples:
        if not isinstance(sample, dict) or not as_bool(sample.get("include")):
            continue
        for field in ("rna_path", "allc_root"):
            if outside_project(output, sample.get(field)):
                data_outside.append("%s:%s" % (sample.get("sample_id", "?"), field))
    if data_outside:
        print("WARNING: these sample inputs live outside the project: %s. Declare them under data.sources in the "
              "intake to have them linked into Data/ instead." % ", ".join(data_outside), file=sys.stderr)
    for warning in reference_warnings:
        print("WARNING: %s" % warning, file=sys.stderr)
    (output / "README.md").write_text(
        bilingual_readme(project_cfg, analysis_cfg["analysis"], scheduler_cfg, samples, env_rows,
                         data_rows, data_outside),
        encoding="utf-8",
    )
    (output / "Report.md").write_text(
        bilingual_report(project_id, provenance, data_rows, data_outside,
                         reference_notes, reference_warnings),
        encoding="utf-8",
    )
    # After the root documents, so a stage document is never adapted before the
    # project it describes exists. The stage READMEs and Reports are copied
    # verbatim from the template; this is the one place that binds them to this
    # project, and the per-stage run logs inside inspect_run.py keep them current
    # from here on.
    stage_labels = [str(sample.get("sample_id")) for sample in samples
                    if isinstance(sample, dict) and as_bool(sample.get("include"))
                    and str(sample.get("sample_id") or "").strip()]
    adapt_stage_docs(output, project_id, project_cfg["project"].get("organism"),
                     project_cfg["project"].get("genome_id"), stage_labels,
                     packaged_notebook(output))
    write_json(output / ".workflow" / "template-lock.json", {
        "skill": "single-cell-multiomics-analysis", "skill_version": provenance["skill_version"],
        "skill_git_commit": provenance["skill_git_commit"], "schema_version": 1,
    })
    print(json.dumps({"project": str(output), "status": "initialized"}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
