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
- 甲基化路线必须显式提供 `chrom_sizes` 和 blacklist；TSS BED 可选。Skill 自动计算 SHA-256，若用户声明期望 checksum，则必须一致；
- DMR 路线需要 `review_status: approved`、非空且唯一的 `cell_id`、无占位标签的 `cell_type`，且至少两个细胞类型达到配置的最小细胞数。

Skill 会先只读发现兼容环境。缺失的 Scanpy/ALLCools、MethSCAn 或 MethylVI 环境会在项目 `.environments/` 下隔离创建、验证并写入 `config/environments.tsv`；不会升级或修改已经存在的共享环境。

### 4. 准备 intake 文件

下面是配对多组学的最小示例。RNA-only 或 ALLC-only 项目可以将另一种输入留空。

`assets/project-template/config/examples/` 另附一份**填好的完整实例** —— `example-ipf-tissue.yaml`（可直接 `init_project.py --intake` 运行）、`example-samples.tsv`，以及记录该项目实测数值的 README。生成的项目不会读取该目录；它的作用是在字段全空之前先看到一个真实填写的样子。其中的 scheduler 段刻意留空：分区、账号、节点名单描述的是你的集群，不是别人的。

```yaml
schema_version: 2
project_id: example_multiome
organism: human

# 使用者数据的当前位置。生成时会把每一项软链接成项目内的 Data/<name>，
# 因此 config/ 里永远不会出现项目外的绝对路径。
data:
  sources:
    - name: Matrix
      path: /shared/lung-tissue/Matrix
    - name: ALLCools
      path: /shared/lung-tissue/ALLCools
    # 也可以直接声明下载来源；生成时不下载，由 tools/link_data.py --execute 取回并校验
    # - name: gse_supplement
    #   url: "https://example.org/gse.tar.gz"
    #   sha256: REPLACE_WITH_SHA256

references:
  genome_id: hg38
  # chrom_sizes 与 blacklist 可以整段省略：技能内置了 GRCh38 的两份公共参考文件，
  # 仅当项目是 human + 有甲基化路线 + genome_id 为 hg38 时自动填入，并在 Report.md 中说明。
  # 要换用自己的文件，就在这里声明（自备 blacklist 需要自己的 blacklist_md5）。
  # chrom_sizes: my_hg38.chrom.sizes
  # blacklist: my_hg38.blacklist.bed
  tss_bed: /reference/hg38_TSS.bed
  tss_bed_sha256: REPLACE_WITH_SHA256

annotation:
  table: Data/annotation/cell_id_cell_type.tsv
  profile: null
  review_status: approved
  cell_id_column: cell_id
  cell_type_column: cell_type

samples:
  - sample_id: sample_A
    condition: control
    batch: batch_1
    include: 1
    rna_path: Data/Matrix/sample_A/filtered_feature_bc_matrix.h5
    rna_format: 10x_h5
    allc_root: Data/ALLCools/sample_A/allc
    allc_glob: "**/*.allc.tsv.gz"
    cell_id_prefix: sample_A
    allc_cell_id_regex: ""
    allc_cell_id_replacement: ""

  - sample_id: sample_B
    condition: treatment
    batch: batch_2
    include: 1
    rna_path: /data/sample_B/filtered_feature_bc_matrix.zip
    rna_format: 10x_zip
    allc_root: /data/sample_B/allc
    allc_glob: "**/*_allc.gz"
    cell_id_prefix: sample_B
    allc_cell_id_regex: ""
    allc_cell_id_replacement: ""

environments:
  - stage: orchestrator
    python: /path/to/analysis-core/bin/python
    executable: /path/to/analysis-core/bin/python
    version_command: /path/to/analysis-core/bin/python --version
    required: 1

  - stage: scanpy_allcools
    python: /path/to/allcools/bin/python
    executable: /path/to/allcools/bin/python
    version_command: /path/to/allcools/bin/python --version
    required: 1

  - stage: methscan
    python: /path/to/methscan/bin/python
    executable: /path/to/methscan/bin/methscan
    version_command: /path/to/methscan/bin/methscan --version
    required: 1

  - stage: methylvi
    python: /path/to/methylvi/bin/python
    executable: /path/to/methylvi/bin/python
    version_command: /path/to/methylvi/bin/python --version
    required: 1

scheduler:
  backend: slurm
  partitions:
    - cpu
  account: null
  allow_nodes: []
  exclude_nodes: []
  max_parallel: 2
  validation_workers: 4
  limited_profiles:
    - methscan_branch
    - dmr
    - feature_builder
    - trainer
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

# 确认每个 Data/ 条目都指向声明的来源（只读；加 --execute 才动手）
"$PYTHON" "$PROJECT/tools/link_data.py" \
  --project "$PROJECT"

# 只读生成环境复用/创建计划
"$PYTHON" "$PROJECT/tools/bootstrap_environments.py" \
  --project "$PROJECT"

# 自动创建缺失的隔离环境并更新 environments.tsv
"$PYTHON" "$PROJECT/tools/bootstrap_environments.py" \
  --project "$PROJECT" \
  --execute

# 环境完成后执行日常快速校验

"$PYTHON" "$PROJECT/tools/validate_project.py" \
  --project "$PROJECT" \
  --mode quick \
  --json-out "$PROJECT/preflight.json"

# 首次正式运行或输入/schema/模板变化后执行全量流式校验
"$PYTHON" "$PROJECT/tools/validate_project.py" \
  --project "$PROJECT" --mode full --workers 8 \
  --json-out "$PROJECT/full-preflight.json"
```

`init_project.py` 会按 `data.sources` 把每份输入软链接成 `Data/<name>`，并把解析后的目标写进 `config/project.yaml`；生成时的状态（已就位 / 待下载 / 本机不可见）同时写进根 `README.md` 与 `Report.md`。生成阶段不下载：`url` 来源由 `tools/link_data.py --execute` 取回，先落到隐藏的临时文件、校验 sha256 通过后才移入 `Data/`，失败时 `Data/` 不留任何残留。`link_data.py` 不会覆盖不是它创建的条目。

`link_data.py` 报告 `dangling` 表示链接本身正确、但本机看不到目标路径——这是挂载问题而非链接问题，`validate_project.py` 会明确拒绝这种条目。若某条路径只在计算节点可见，就在执行上下文里校验。

bootstrap 优先复用通过导入测试的已有环境；缺失 profile 默认创建到 `PROJECT/.environments/`。计划和执行证据保存在 `.workflow/environment-bootstrap/`。如需使用其他可写位置，可传入 `--prefix-root`。

quick 模式逐文件检查索引、格式和有界记录；full 模式并行流式检查全部 ALLC 的列、context、排序和 `mc <= cov`。正式提交要求与当前输入签名一致的 full 证据；两种校验都不会修改输入数据或共享环境。

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

`plan_workflow.py` 只生成 DAG，不提交作业。`--routes auto` 根据输入启用可运行路线，也可以显式传入逗号分隔的路线名。显式选择子路线时，只校验该 DAG 实际使用的输入模态、参考文件和环境；例如 `--routes scanpy` 不扫描 ALLC，也不要求 MethSCAn/MethylVI 环境或甲基化参考文件。ALLCools 自动补齐到 MethSCAn filter 的基础链，但不会因此额外运行 smooth/VMR。

Scanpy、MethSCAn、ALLCools 和 MethylVI 均有内置执行适配器，并自动使用 `environments.tsv` 中对应阶段的 Python。`analysis.task_commands` 仅用于专家覆盖；覆盖命令成功后必须写出计划声明的 `<task_dir>/task_outputs.json`：

```json
{
  "artifacts": ["/absolute/path/to/result"]
}
```

`artifacts` 必须是非空列表，且每个路径都必须存在。任务状态会分别保存子进程退出码和证据校验后的包装任务退出码。

### 7. 正式提交

确认 preflight、DAG、task commands 和 dry-run 后，移除 `--dry-run`：

```bash
"$PYTHON" "$PROJECT/tools/submit_workflow.py" \
  --project "$PROJECT" \
  --run-id run_001
```

每个 Slurm 任务在执行 `sbatch` 前都会重新读取 `sinfo`、`scontrol show nodes`、全局 `squeue` 和当前用户 `squeue`。资源报告分别保存 Slurm 可调度 CPU/内存和操作系统观测到的 `CPULoad`/`FreeMem`。

当资源低于正确性所需 floor、分区不存在或 Slurm 查询失败时，不会强行提交。默认按 partition 调度，不固定节点。local 后端使用同一个 DAG，但受到 `max_threads`、`max_memory` 和并发设置限制。

所有正式结果写入 `PROJECT/Results/runs/<run_id>/`。同一 run 中断后可再次执行 submit：已完成或仍在调度器中的任务不会重复提交，失败或产物缺失的任务会按依赖恢复。

### 8. 监控和更新报告

```bash
"$PYTHON" "$PROJECT/tools/inspect_run.py" \
  --project "$PROJECT" \
  --run-id run_001

# 需要时（清理 run 目录后、从旧版本迁移后、区域被手工改过后）重建整个运行记录：
"$PYTHON" "$PROJECT/tools/update_report.py" \
  --project "$PROJECT"
```

`inspect_run.py` 在给出运行结论的同时，会把本次 run 的精简记录写入根 `Report.md` 的运行记录区域，因此**多次迭代会累积**而不是互相覆盖；重跑同一 run 会就地更新它自己的记录。同一条记录也会写入本次 run 触及的每个阶段 Report（`Scripts/Scanpy/Report.md` 等）各自的运行记录区域，内容是同一份 summary 按该阶段的任务过滤后的视图 —— 根日志始终是完整记录。这条写入是尽力而为：失败只打印 stderr 警告并附上恢复命令，不改变 run 的退出码。`update_report.py` 是显式路径，写不进去就非零退出；不带 `--run-id` 时从 `.workflow/runs/*/run_summary.json` 全量重建根日志与全部阶段日志（`--run-id` 只用于要求该 run 必须存在）。

根 `README.md` 由 `init_project.py` 依据 intake 与配置**写入一次后冻结**，不反映进度或结果；`Scripts/<Module>/` 下的模块 README 与 Report 正文同样由生成时**一次写入**——`init_project.py` 会在每个阶段文档顶部写入一段标明本项目（项目名、物种、基因组、实际样本、notebook 入口）的上下文块，并把 `<notebook>` 占位符解析为本项目实际入口。阶段 Report 正文中出现的示例样本名与实测数值属于参考项目，由该上下文块明确标注；运行只写入其运行记录区域，不改写正文。

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

对应分析结果位于 `PROJECT/Results/runs/<run_id>/`，不同 run 不共享可写结果目录。

不能只根据作业已经提交或 Slurm 显示 `COMPLETED` 判断分析完成。任务输出、签名、机器可读 summary 和完成标记必须一致。

### 9. 审核并记录细胞类型

Scanpy 给出的 cell type 是**候选**；依赖 cell type 的 DMR 路线只有拿到人工审核过的 profile 才会被规划。审核走两步，第一步生成逐聚类的审核表（`cell_type` 一列已用该聚类自身占多数的标签预填），改完再记录：

```bash
# 1. 生成审核表（一个聚类一行）
"$PYTHON" "$PROJECT/tools/record_annotation_review.py" \
  --project "$PROJECT" --run-id run_001 --worksheet review.tsv

# 2. 修改 cell_type 一列后记录
"$PYTHON" "$PROJECT/tools/record_annotation_review.py" \
  --project "$PROJECT" --run-id run_001 --mapping review.tsv
```

工具会把该 run 的 `analysis_signature` 与聚类集合一起写进 profile（默认 `config/annotation.yaml`，或 `annotation.profile` 指定的路径），因此改过 `leiden.resolution` 或换了数据后旧 profile 自动失效——**重新审核**，不要把映射搬过去。占位标签（`NA`、`Unassigned`、`requires_review` 等）、覆盖不全的审核表、以及会被覆盖的既有审核表都会被拒绝。记录后还需在 `config/project.yaml` 里设置 `annotation.profile`，并把 `annotation.review_status` 改为 `approved`；两者齐备 DMR 路线才会被规划。

**注释在产出它的那一轮里完成，不停在 `Unassigned`。** candidate 运行的作用是产出证据（逐 cluster 的排名 marker、marker 图、QC），判定标签是同一轮的下一步：读该轮自己的 marker 证据定出每个 cluster 的 cell type，写进 `markers.dotplot_markers` 与 `markers.cell_type_order`，再用上面的工具记录，然后以 `baseline` 重跑。skill 不内置 marker 面板、也不提供对照表——沿用别的项目的标签正是这套流程要避免的移植。标签定下来后，交付的点图是 baseline 的 `annotation_marker_dotplot.png`：**每行一个细胞类型、样本合并在行内**；candidate 那张纵轴为 cluster 的 `candidate_cluster_marker_dotplot.png` 是判定时读的证据，会被它取代。最后请用户确认注释，说明哪些 cluster 已确定、哪些不确定。详见 [Scanpy 与 Harmony](references/scanpy.md)。

### 10. 测试

```bash
PYTHON_EXE=/path/to/python-3.9-or-newer bash tests/run_tests.sh

/path/to/python-3.9-or-newer \
  "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" \
  "$HOME/single-cell-multiomics-analysis"
```

仓库测试覆盖 RNA-only、ALLC-only、配对多组学、环境 bootstrap、输入数据链接与下载校验、内置参考默认值、完成证据、路线闭包、常见 ALLC 命名、错误 checksum、Slurm 资源和 local backend。

### 11. 进一步文档

- [配置与 intake](references/configuration.md)
- [环境发现与自动配置](references/environments.md)
- [Scanpy 路线](references/scanpy.md)
- [MethSCAn 路线](references/methscan.md)
- [ALLCools 和 MethylVI 路线](references/methylvi.md)
- [调度与执行](references/execution.md)
- [失败恢复](references/troubleshooting.md)
- [文档与运行证据](references/documentation.md)
- [模板维护与迁移清单](references/template-maintenance.md)、[迁移记录](references/migration-inventory.tsv)

模板内每个模块另有一份完整的中英双语说明与证据台账，随项目一起生成：

- 项目根 `README.md`（生成时一次写入、此后冻结的稳定契约）与 `Report.md`（按 run 累积运行记录，区域格式见 [documentation](references/documentation.md)）
- `Scripts/Scanpy/README.md`、`Scripts/Scanpy/Report.md`
- `Scripts/Methscan/README.md`、`Scripts/Methscan/Report.md`
- `Scripts/Methylvi/README.md`、`Scripts/Methylvi/Report.md`，以及 `allcools/`、`vmr/`、`vmr_dmr/`、`shared/` 各自的 README
- `Scripts/Environment/README.md`、`Scripts/Environment/Report.md`
- `config/examples/README.md` —— 示例项目的实测数值，以及哪些值必须重新推导而不可照搬
- `Supplementary/README.md` —— 内置参考文件的来源、摘要，以及如何替换为自己的文件

`Scripts/Scanpy/Notebooks/` 下有两个 notebook：`scanpy_workflow.ipynb` 是 run DAG 实际执行的入口，`example_ipf_scanpy.ipynb` 是随包附带、**不参与执行**的完整实例，用于展示一次真实的端到端探索（其中的人工注释映射只对那一组 cluster 有效，不可移植）。

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

The intake's `data` block declares where the user's data already lives; generation links each entry into the project's own `Data/`, and the sample paths are then written relative to the project root (`Data/Matrix/...`). Nothing under `config/` points outside the project. `references.chrom_sizes` and `references.blacklist` may be omitted entirely: the skill ships public GRCh38 files under `Supplementary/` and fills in whichever key is blank for a human project with a methylation route — and says so in the generated `Report.md`. A project on another organism or genome keeps the preflight failure that names the reference it must declare, rather than silently inheriting the wrong coordinates.

```bash
SKILL_ROOT="$HOME/single-cell-multiomics-analysis"
PYTHON=/path/to/python-3.9-or-newer
PROJECT=/work/example_multiome

"$PYTHON" "$SKILL_ROOT/scripts/init_project.py" \
  --intake intake.yaml \
  --output "$PROJECT"

"$PYTHON" "$PROJECT/tools/link_data.py" \
  --project "$PROJECT" \
  --execute

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

Generation links every declared `data.sources` entry into the project's own `Data/<name>` and writes the resolved target into `config/project.yaml`; the state at generation time (in place, pending download, or not visible from this host) also goes into the root `README.md` and `Report.md`. Generation never downloads: a `url` source is fetched by `tools/link_data.py --execute`, which streams it to a hidden temporary file, verifies the declared sha256, and only then moves it into `Data/`. A failed fetch leaves nothing behind. `link_data.py` never overwrites an entry it did not create.

A `dangling` state from `link_data.py` means the link is exactly right but its target is unreachable from this host — a mount question, not a link question. `validate_project.py` refuses such an entry by name, so check data reachability in the execution context.

The bootstrapper reuses compatible environments read-only and creates missing profiles under `PROJECT/.environments/` from bundled versioned specs. It verifies imports or executables before updating `config/environments.tsv`; it never upgrades a discovered shared environment.

Review the preflight report, DAG, and resource plans. Built-in adapters cover Scanpy, MethSCAn, ALLCools, and MethylVI and select the declared stage environment. Explicit route subsets validate only the input modalities, references, and environments used by their closed DAG. `analysis.task_commands` is an expert override, not a required setup step. A successful override must create `<task_dir>/task_outputs.json` with a non-empty `artifacts` array of existing absolute paths; child-process and wrapper-task return codes are recorded separately.

To execute an approved run, repeat `submit_workflow.py` without `--dry-run`. Production submission requires matching full-validation evidence and creates it when absent. Each Slurm task receives a fresh resource query immediately before submission. Local execution uses the same DAG under configured limits. Results are isolated under `Results/runs/<run_id>/`, and re-submission resumes incomplete work.

`assets/project-template/config/examples/` carries a completed intake (`example-ipf-tissue.yaml`, runnable directly with `init_project.py --intake`), the same rows as a generated manifest, and a README recording one study's observed magnitudes and the values that must be re-derived rather than copied. Generated projects never read it. Its scheduler block is deliberately empty: partitions, accounts, and node lists describe your cluster, not someone else's.

The template also ships standalone `run_*.sbatch` and `submit_*.sh` wrappers that run one stage without a plan, for exploring a stage, smoke testing a change, or recovering a single stage of a completed run. They resolve resources from the same named scheduler profiles the DAG uses and carry no site directive, but they produce no `task_outputs.json` and are never a run's completion evidence. The MethylVI route runners additionally require an explicit `SCMO_STANDALONE_ACK=1`, which those wrappers set; a direct `run.sh` call is never granted it implicitly.

### Further documentation

- Skill-level references: [configuration](references/configuration.md), [environments](references/environments.md), [Scanpy](references/scanpy.md), [MethSCAn](references/methscan.md), [ALLCools and MethylVI](references/methylvi.md), [execution](references/execution.md), [troubleshooting](references/troubleshooting.md), [documentation](references/documentation.md), [template maintenance](references/template-maintenance.md).
- The project root carries `README.md` — the stable contract, written once at generation and frozen — and `Report.md`, which accumulates one concise record per run in a tool-owned region (format in [documentation](references/documentation.md)).
- Every packaged module also carries its own bilingual README and evidence ledger, generated with the project: `Scripts/Scanpy/`, `Scripts/Methscan/`, `Scripts/Methylvi/` (plus `allcools/`, `vmr/`, `vmr_dmr/`, `shared/`), and `Scripts/Environment/`. Generation writes each one once, adding a context block that names this project, its organism, genome, actual samples and notebook entry; the prose below it — including the reference project's example sample names and worked numbers — is labelled as the reference project's and is never rewritten. The stage Report additionally carries its own run-log region, refreshed per run with that stage's tasks only, so a reader of `Scripts/Scanpy/Report.md` does not have to read the MethylVI runs around it.
- `Supplementary/README.md` records the provenance and digests of the bundled reference files, and how to override them.
- `Scripts/Scanpy/Notebooks/scanpy_workflow.ipynb` is the notebook the run DAG executes. `example_ipf_scanpy.ipynb` beside it is a complete worked example that is never executed; its reviewed cluster-to-cell-type mapping is valid only for that example's exact cluster set and must not be transplanted.

### Inspect completion

```bash
"$PYTHON" "$PROJECT/tools/inspect_run.py" \
  --project "$PROJECT" \
  --run-id run_001

# When needed — after pruning run directories, migrating an older project, or hand edits:
"$PYTHON" "$PROJECT/tools/update_report.py" \
  --project "$PROJECT"
```

`inspect_run.py` records the run's verdict and, at the same time, refreshes that run's concise entry in the root `Report.md` run log, so successive iterations **accumulate** instead of overwriting each other; re-inspecting a run updates its own entry in place. The same entry also goes into the run-log region of every stage Report the run touched, filtered to that stage's tasks — a per-stage view of one corpus of summaries, not a second source of truth. That write is best effort: a failure prints a stderr warning with the recovery command and leaves the run's exit code untouched. `update_report.py` is the loud path and exits non-zero when it cannot write; without `--run-id` it rebuilds the root log and every stage log from `.workflow/runs/*/run_summary.json` (`--run-id` only requires that one run to exist).

The root `README.md` is written once by `init_project.py` from the intake and configuration and is **frozen** afterwards: it reports no progress and no results. The module READMEs and the prose of the module Reports under `Scripts/<Module>/` are written once at generation too, each carrying a generated context block naming this project, its organism, genome, actual samples and notebook entry. The example sample names and measured numbers in that prose belong to the reference project and are labelled as such; rewriting them into statements about this project would fabricate evidence. Only the stage Reports' run-log regions are updated by a run.

A run is complete only when validated outputs, signatures, machine-readable summaries, and completion markers agree. Slurm submission or a `COMPLETED` scheduler state alone is not sufficient.

### Review and record cell types

Scanpy's cell types are **proposals**. The cell-type DMR routes are planned only against a profile a human has reviewed, and the review is two passes: the first writes a per-cluster worksheet whose `cell_type` column is pre-filled with that cluster's own majority label, and the second records what you changed.

```bash
# 1. a worksheet, one row per cluster
"$PYTHON" "$PROJECT/tools/record_annotation_review.py" \
  --project "$PROJECT" --run-id run_001 --worksheet review.tsv

# 2. after editing the cell_type column
"$PYTHON" "$PROJECT/tools/record_annotation_review.py" \
  --project "$PROJECT" --run-id run_001 --mapping review.tsv
```

The tool writes the run's `analysis_signature` and cluster set into the profile alongside the labels (to `config/annotation.yaml`, or to `annotation.profile` if it names a path), so a profile recorded for one set of parameters cannot be reused after `leiden.resolution` changes or the data is re-run — review again instead of transplanting the mapping. Placeholder labels (`NA`, `Unassigned`, `requires_review`, and similar), a worksheet that does not cover the run's clusters exactly, and an existing worksheet it would overwrite are all refused. Recording is not enough on its own: set `annotation.profile` in `config/project.yaml` and set `annotation.review_status` to `approved`, and only then do the DMR routes become plannable.

**The first pass annotates; it does not stop at `Unassigned`.** A candidate run exists to produce the evidence — per-cluster ranked markers, the marker figure, QC — and deciding the labels from it is the next step of the same pass: read that run's own marker evidence, settle a cell type per cluster, write the labels into `markers.dotplot_markers` and `markers.cell_type_order`, record them with the tool above, and re-run as `baseline`. No marker panel ships with the skill and there is no crosswalk to consult, because carrying another project's labels over is the transplant this workflow exists to prevent. Once the labels are recorded, the dotplot that ships is the baseline's `annotation_marker_dotplot.png` — **one row per cell type, samples pooled inside each row**; the candidate's `candidate_cluster_marker_dotplot.png`, whose rows are Leiden clusters, was the evidence read while deciding and is superseded by it. Close by asking the user to confirm the annotation, naming the clusters you consider settled and the ones you do not. See [Scanpy and Harmony](references/scanpy.md).

See the linked references above for configuration, route-specific behavior, resource selection, troubleshooting, and evidence contracts.
