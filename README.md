# Single-cell Multiomics Analysis

A portable Codex skill for configuration-driven single-cell RNA and DNA-methylation analysis with Scanpy/Harmony, MethSCAn, ALLCools, and MethylVI.

一个由配置驱动的通用 Codex Skill，用于执行 Scanpy/Harmony 单细胞 RNA、MethSCAn VMR/DMR、ALLCools 特征构建和 MethylVI 分析。

## 中文使用说明

### 1. 功能范围

支持以下输入组合：

| 输入 | 可用路线 |
|---|---|
| 仅 RNA | Scanpy/Harmony、QC、聚类、marker 审核表和 UMAP |
| 仅逐细胞 ALLC | MethSCAn VMR、ALLCools、无监督 MethylVI |
| RNA + ALLC | RNA 注释、VMR、分样本/合并样本 DMR、ALLCools、MethylVI |
| ALLC + 有效 `cell_type` | pairwise/pooled DMR、三种热图、VMR + 全部 unique pooled-DMR MethylVI |

RNA 输入支持标准 10x filtered matrix 目录、ZIP 和 10x H5。甲基化输入支持逐细胞 `.allc.tsv.gz`，也支持通过显式 glob 选择 ALLCools 风格的 `*_allc.gz`；每个压缩 ALLC 必须具有相邻的 `.tbi`。

本项目不处理原始 FASTQ 或 BAM，也不内置任何物种、基因组、样本名、细胞数、固定节点或用户目录。

### 2. 安装

```bash
git clone git@github.com:Zhixiong26/single-cell-multiomics-analysis.git \
  "$HOME/single-cell-multiomics-analysis"

mkdir -p "$HOME/.codex/skills"
ln -s "$HOME/single-cell-multiomics-analysis" \
  "$HOME/.codex/skills/single-cell-multiomics-analysis"
```

如果目标软链接已经存在，请先确认它指向的仓库，不要直接覆盖。安装后开启新的 Codex 会话；可通过以下提示词调用：

```text
请使用 $single-cell-multiomics-analysis，为我提供的 10x RNA 和逐细胞 ALLC
创建一个新项目。先生成 intake、执行完整 preflight 和 dry-run，不要立即提交正式任务。
```

生成后的项目包含自己的脚本和运行工具，即使之后没有安装本 Skill，也可独立运行。

### 3. 前置条件

- orchestrator Python 3.9 或更高版本；
- 对输入数据和参考文件具有只读权限；
- 系统中具有 `conda` 或 `mamba`；已有兼容分析环境可以复用，没有时由 Skill 创建；
- Slurm 模式下可以运行 `sinfo`、`scontrol`、`squeue`、`sbatch` 和 `sacct`；
- genome manifest 显式提供 `chrom_sizes`、blacklist 和可选 TSS BED，推荐同时提供 SHA-256；
- DMR 路线需要非空、唯一的 `cell_id` 和有效的 `cell_type`，且至少两个细胞类型达到配置的最小细胞数。

Skill 会先只读发现兼容环境。缺失的 Scanpy/ALLCools、MethSCAn 或 MethylVI 环境会在项目 `.environments/` 下隔离创建、验证并写入 `config/environments.tsv`；不会升级或修改已经存在的共享环境。

### 4. 准备 intake 文件

下面是配对多组学的最小示例。RNA-only 或 ALLC-only 项目可以将另一种输入留空。

```yaml
schema_version: 1
project_id: example_multiome
organism: human

references:
  genome_id: hg38
  chrom_sizes: /reference/hg38.chrom.sizes
  chrom_sizes_sha256: REPLACE_WITH_SHA256
  blacklist: /reference/hg38.blacklist.bed
  blacklist_sha256: REPLACE_WITH_SHA256
  tss_bed: /reference/hg38_TSS.bed
  tss_bed_sha256: REPLACE_WITH_SHA256

annotation:
  path: /data/annotation/cell_id_cell_type.tsv
  cell_id_column: cell_id
  cell_type_column: cell_type

samples:
  - sample_id: sample_A
    condition: control
    batch: batch_1
    include: 1
    rna_path: /data/sample_A/filtered_feature_bc_matrix.h5
    rna_format: 10x_h5
    allc_root: /data/sample_A/allc
    allc_glob: "**/*.allc.tsv.gz"
    cell_id_prefix: sample_A

  - sample_id: sample_B
    condition: treatment
    batch: batch_2
    include: 1
    rna_path: /data/sample_B/filtered_feature_bc_matrix.zip
    rna_format: 10x_zip
    allc_root: /data/sample_B/allc
    allc_glob: "**/*_allc.gz"
    cell_id_prefix: sample_B

environments:
  - stage: orchestrator
    python: /envs/allcools/bin/python
    executable: /envs/allcools/bin/python
    version_command: /envs/allcools/bin/python --version
    required: 1

  - stage: scanpy_allcools
    python: /envs/allcools/bin/python
    executable: /envs/allcools/bin/python
    version_command: /envs/allcools/bin/python --version
    required: 1

  - stage: methscan
    python: /envs/methscan/bin/python
    executable: /envs/methscan/bin/methscan
    version_command: /envs/methscan/bin/methscan --version
    required: 1

  - stage: methylvi
    python: /envs/methylvi/bin/python
    executable: /envs/methylvi/bin/python
    version_command: /envs/methylvi/bin/python --version
    required: 1

scheduler:
  backend: slurm
  partitions:
    - cpu
  account: null
  allow_nodes: []
  exclude_nodes: []
  max_parallel: 2
  memory_headroom_mb: 4096

analysis:
  methscan:
    min_cells: 6
    dmr:
      raw_p: 0.01
      min_abs_diff: 0.25
      top_per_cell_type: 200
  methylvi:
    feature_targets:
      - 10000
      - 30000
    validation_fraction: 0.1
```

`top_per_cell_type` 只用于 Top-N hypo-DMR 热图。VMR+DMR MethylVI 路线使用满足阈值的全部 unique pooled-DMR，两者不会混用。详细字段见 [configuration.md](references/configuration.md)。

如果使用者没有现成环境，可以省略 `environments` 或保留待创建的目标路径，后续 bootstrap 会自动补齐。

### 5. 生成项目并配置环境

输出目录必须不存在或为空：

```bash
SKILL_ROOT="$HOME/single-cell-multiomics-analysis"
PYTHON=/path/to/python-3.9-or-newer
PROJECT=/work/example_multiome

"$PYTHON" "$SKILL_ROOT/scripts/init_project.py" \
  --intake intake.yaml \
  --output "$PROJECT"

# 先只读生成环境复用/创建计划
"$PYTHON" "$PROJECT/tools/bootstrap_environments.py" \
  --project "$PROJECT"

# 自动创建缺失的隔离环境并更新 environments.tsv
"$PYTHON" "$PROJECT/tools/bootstrap_environments.py" \
  --project "$PROJECT" \
  --execute

# 环境完成后执行完整校验

"$PYTHON" "$PROJECT/tools/validate_project.py" \
  --project "$PROJECT" \
  --json-out "$PROJECT/preflight.json"
```

bootstrap 优先复用通过导入测试的已有环境；缺失 profile 默认创建到 `PROJECT/.environments/`。计划和执行证据保存在 `.workflow/environment-bootstrap/`。如需使用其他可写位置，可传入 `--prefix-root`。

校验内容包括配置 schema、文件和 `.tbi`、ALLC 列和 CG context、`mc <= cov`、cell ID、cell type、参考文件 checksum、环境版本和路线依赖。校验不会修改输入数据或共享环境。

### 6. 规划 DAG 和 dry-run

```bash
"$PYTHON" "$PROJECT/tools/plan_workflow.py" \
  --project "$PROJECT" \
  --routes auto \
  --run-id run_001

"$PYTHON" "$PROJECT/tools/submit_workflow.py" \
  --project "$PROJECT" \
  --run-id run_001 \
  --dry-run
```

`plan_workflow.py` 只生成 DAG，不提交作业。`--routes auto` 根据输入启用可运行路线，也可以显式传入逗号分隔的路线名。

Scanpy 已包含默认执行适配器。正式运行 MethSCAn、ALLCools 或 MethylVI 前，需要由 Codex 根据用户的环境和输入在 `config/analysis.yaml` 中生成相应的 `analysis.task_commands`。缺少命令时，生产提交会在任何 `sbatch` 之前被阻止；dry-run 仍可用于审核 DAG 和资源方案。

### 7. 正式提交

确认 preflight、DAG、task commands 和 dry-run 后，移除 `--dry-run`：

```bash
"$PYTHON" "$PROJECT/tools/submit_workflow.py" \
  --project "$PROJECT" \
  --run-id run_001
```

每个 Slurm 任务在执行 `sbatch` 前都会重新读取 `sinfo`、`scontrol show nodes`、全局 `squeue` 和当前用户 `squeue`。资源报告分别保存 Slurm 可调度 CPU/内存和操作系统观测到的 `CPULoad`/`FreeMem`。

当资源低于正确性所需 floor、分区不存在或 Slurm 查询失败时，不会强行提交。默认按 partition 调度，不固定节点。local 后端使用同一个 DAG，但受到 `max_threads`、`max_memory` 和并发设置限制。

### 8. 监控和更新报告

```bash
"$PYTHON" "$PROJECT/tools/inspect_run.py" \
  --project "$PROJECT" \
  --run-id run_001

"$PYTHON" "$PROJECT/tools/update_report.py" \
  --project "$PROJECT" \
  --run-id run_001
```

每次运行的证据位于：

```text
PROJECT/.workflow/runs/<run_id>/
├── plan.json
├── preflight.json
├── submissions.json
├── resource_snapshots/
├── resource_usage/
├── tasks/
├── run_summary.json
└── workflow.COMPLETE
```

不能只根据作业已经提交或 Slurm 显示 `COMPLETED` 判断分析完成。任务输出、签名、机器可读 summary 和完成标记必须一致。

### 9. 测试

```bash
/envs/allcools/bin/python -m unittest discover \
  -s tests \
  -p 'test_*.py' \
  -v

/envs/allcools/bin/python \
  "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" \
  "$HOME/single-cell-multiomics-analysis"
```

仓库测试覆盖 RNA-only、ALLC-only、配对多组学、环境 bootstrap 计划、重复 cell ID、错误 checksum、缺失环境、Slurm 节点状态、资源不足和 local fallback。

### 10. 进一步文档

- [配置与 intake](references/configuration.md)
- [环境发现与自动配置](references/environments.md)
- [Scanpy 路线](references/scanpy.md)
- [MethSCAn 路线](references/methscan.md)
- [ALLCools 和 MethylVI 路线](references/methylvi.md)
- [调度与执行](references/execution.md)
- [失败恢复](references/troubleshooting.md)
- [文档与运行证据](references/documentation.md)

## English guide

### Scope

This skill generates auditable projects from standard 10x RNA matrices and indexed per-cell ALLC files. It supports RNA-only, ALLC-only, and paired multiomics projects. Available routes include Scanpy/Harmony, MethSCAn VMR, per-sample and pooled DMR, hypo-DMR heatmaps, ALLCools feature construction, baseline MethylVI, and VMR plus all unique pooled-DMR MethylVI.

Raw FASTQ and BAM processing are out of scope. No genome, sample label, cell count, filesystem path, Slurm node, or cluster annotation is assumed.

### Install

```bash
git clone git@github.com:Zhixiong26/single-cell-multiomics-analysis.git \
  "$HOME/single-cell-multiomics-analysis"

mkdir -p "$HOME/.codex/skills"
ln -s "$HOME/single-cell-multiomics-analysis" \
  "$HOME/.codex/skills/single-cell-multiomics-analysis"
```

Start a new Codex session and invoke it explicitly when desired:

```text
Use $single-cell-multiomics-analysis to generate a project for my 10x RNA and
per-cell ALLC inputs. Create the intake, run preflight and a dry-run, and do not
submit production jobs yet.
```

### Run

Prepare an intake file using the schema shown above, then run:

```bash
SKILL_ROOT="$HOME/single-cell-multiomics-analysis"
PYTHON=/path/to/python-3.9-or-newer
PROJECT=/work/example_multiome

"$PYTHON" "$SKILL_ROOT/scripts/init_project.py" \
  --intake intake.yaml \
  --output "$PROJECT"

"$PYTHON" "$PROJECT/tools/bootstrap_environments.py" \
  --project "$PROJECT" \
  --execute

"$PYTHON" "$PROJECT/tools/validate_project.py" \
  --project "$PROJECT" \
  --json-out "$PROJECT/preflight.json"

"$PYTHON" "$PROJECT/tools/plan_workflow.py" \
  --project "$PROJECT" \
  --routes auto \
  --run-id run_001

"$PYTHON" "$PROJECT/tools/submit_workflow.py" \
  --project "$PROJECT" \
  --run-id run_001 \
  --dry-run
```

The bootstrapper reuses compatible environments read-only and creates missing profiles under `PROJECT/.environments/` from bundled versioned specs. It verifies imports or executables before updating `config/environments.tsv`; it never upgrades a discovered shared environment.

Review the preflight report, DAG, task commands, and resource plans. Scanpy has a default execution adapter. MethSCAn, ALLCools, and MethylVI production tasks require commands adapted to the declared inputs and environments in `config/analysis.yaml`; submission is blocked before the first `sbatch` if any planned command is missing.

To execute an approved run, repeat `submit_workflow.py` without `--dry-run`. Each Slurm task receives a fresh resource query immediately before submission. Local execution uses the same DAG under configured thread, memory, and concurrency limits.

### Inspect completion

```bash
"$PYTHON" "$PROJECT/tools/inspect_run.py" \
  --project "$PROJECT" \
  --run-id run_001

"$PYTHON" "$PROJECT/tools/update_report.py" \
  --project "$PROJECT" \
  --run-id run_001
```

A run is complete only when validated outputs, signatures, machine-readable summaries, and completion markers agree. Slurm submission or a `COMPLETED` scheduler state alone is not sufficient.

See the linked references above for configuration, route-specific behavior, resource selection, troubleshooting, and evidence contracts.
