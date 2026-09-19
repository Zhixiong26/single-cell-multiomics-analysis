# VMR + all unique pooled-DMR MethylVI / VMR + 全部 unique pooled-DMR MethylVI

This workflow adds pooled cell-type DMR information to each existing
MethSCAn-VMR MethylVI input and retrains six models without modifying the
completed baseline results.

本流程为每个已有的 MethSCAn-VMR MethylVI 输入加入 pooled cell-type DMR 信息，
并重新训练 6 个模型，不修改已完成的基线结果。

## Fixed analysis definition / 固定分析定义

- Source: every pooled cell-type comparison produced by the completed pooled-DMR
  job under `Results/Methscan/<run-name>/pooled/07_methdiff` (reference project:
  all 91 comparisons among 14 cell types).
- Retain every row with raw `p < 0.01` and `abs(meth_A - meth_B) >= 0.25`.
- Define the hypo cell type from MethSCAn's `low_group_label`.
- Remove exact duplicates within each hypo cell type and merge all overlapping
  retained DMR intervals. There is no Top-N/Top200 truncation.
- Apply the same blacklist rule used by the VMR workflow — remove source
  intervals with at least 20% blacklist overlap (`blacklist_fraction=0.2` on the
  configured reference blacklist).
- Preserve each existing Top10k/Top30k VMR input. Append only merged DMRs that
  have zero genomic overlap with that model's selected VMRs, preventing the same
  CpG from being represented twice.
- Train `{0.01,0.02,0.05} x {Top10k,Top30k}` using integer `mc/cov`, the same
  `$SCMO_EXPECTED_CELLS` cells, `sample_id` batch key, seed 0, and early
  stopping.

- 来源：已完成 pooled-DMR 作业在 `Results/Methscan/<run-name>/pooled/07_methdiff`
  下的全部 pooled cell-type 比较（参考项目：14 个 cell type 之间的全部 91 组
  比较）。
- 保留所有满足 raw `p < 0.01` 与 `abs(meth_A - meth_B) >= 0.25` 的行。
- 依据 MethSCAn 的 `low_group_label` 定义 hypo cell type。
- 在每个 hypo cell type 内去除完全重复的行，并合并所有重叠的保留 DMR 区间。
  不做 Top-N/Top200 截断。
- 应用与 VMR 流程相同的 blacklist 规则——移除与 blacklist 重叠至少 20% 的源
  区间（在配置的参考 blacklist 上 `blacklist_fraction=0.2`）。
- 保留每个已有的 Top10k/Top30k VMR 输入。只追加与该模型所选 VMR 零基因组重叠
  的合并 DMR，避免同一 CpG 被重复表示。
- 使用整数 `mc/cov`、同一批 `$SCMO_EXPECTED_CELLS` 个细胞、`sample_id` batch
  key、seed 0 与 early stopping 训练 `{0.01,0.02,0.05} x {Top10k,Top30k}`。

The completed preparation currently contains 1,104,139 exact-unique hypo-DMRs
and 321,910 merged non-overlapping DMR intervals. The merged regions cover
967,355,000 bp and the maximum interval length is 58 kb.

已完成的准备步骤目前包含 1,104,139 个完全唯一的 hypo-DMR 与 321,910 个合并后的
互不重叠 DMR 区间。合并区域覆盖 967,355,000 bp，最大区间长度为 58 kb。

## Submission / 提交

Run from the project root (`$SCMO_PROJECT_ROOT`):

在项目根目录（`$SCMO_PROJECT_ROOT`）下运行：

```bash
bash Scripts/Methylvi/vmr_dmr/submit_vmr_dmr_pipeline.sh
```

The submitter reuses the completed DMR preparation, builds the pooled-DMR count
matrix once, derives six model-specific non-overlap unions, trains the six
models, generates the standard cell-type/sample/condition/Leiden UMAPs, and
validates all outputs before writing `workflow.COMPLETE`.

提交脚本复用已完成的 DMR 准备结果，只构建一次 pooled-DMR 计数矩阵，派生 6 个
模型专属的非重叠并集，训练 6 个模型，生成标准的
cell-type/sample/condition/Leiden UMAP，并在写出 `workflow.COMPLETE` 之前校验
全部输出。

All jobs run on a single node (node selection comes from the scheduler profile,
not from this document). Count construction uses 40 CPU/48G; input joins use
4 CPU/40G; training uses 32 CPU/48G. The six join/train pairs are serialized to
keep only one memory-intensive task active on the node at a time.

所有作业都运行在同一个节点上（节点选择由 scheduler profile 决定，而不是由本文档
决定）。计数构建使用 40 CPU/48G；输入 join 使用 4 CPU/40G；训练使用
32 CPU/48G。6 组 join/train 串行执行，以保证同一时刻该节点上只有一个高内存
任务在运行。

## Current production run / 当前正式运行

Status date / 状态日期: `<date>`. DMR preparation is reused; the count job
`<job-id>` is running. Serialized join/train jobs are `<job-id>`–`<job-id>`, and
the final validator `<job-id>` waits for all six models. The workflow remains
`in progress` until `run_summary.json` and `workflow.COMPLETE` are created.

状态日期：`<date>`。DMR 准备结果被复用；计数作业 `<job-id>` 正在运行。串行化的
join/train 作业为 `<job-id>`–`<job-id>`，最终校验作业 `<job-id>` 等待全部 6 个
模型。在 `run_summary.json` 与 `workflow.COMPLETE` 创建之前，流程保持
`in progress`。

Results are written only below the VMR+DMR run root:

结果只写入 VMR+DMR run 根目录之下：

```text
Results/MethylVI/<run-name>/
```

This root is separate from the completed baseline root, which is never
overwritten.

该根目录与已完成的基线根目录相互独立，基线永不被覆盖。
