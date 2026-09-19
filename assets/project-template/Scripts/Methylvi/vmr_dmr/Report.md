# VMR + all unique pooled-DMR MethylVI report / VMR + 全部 unique pooled-DMR MethylVI 报告

## Input definition / 输入定义

- Cells: `$SCMO_EXPECTED_CELLS` MethSCAn-filtered cells (reference project:
  6,264 = CYL 2,919 + ZCP 3,345; sample names are examples only).
- Pooled DMR source: 91/91 completed comparisons from the MethSCAn pooled job
  `<job-id>` (run `Results/Methscan/<run-name>`).
- DMR filter: raw `p < 0.01` and `abs(methdiff) >= 0.25`.
- No Top200 or other Top-N truncation.
- Exact deduplication key: hypo cell type + chromosome + start + end.
- Blacklist: remove source intervals with at least 20% overlap.
- DMR merge: merge all overlapping retained intervals.
- VMR/DMR overlap: retain selected VMRs and append only zero-overlap DMRs.

- 细胞：`$SCMO_EXPECTED_CELLS` 个 MethSCAn filter 后的细胞（参考项目：6,264，即
  CYL 2,919 + ZCP 3,345；样本名仅为示例）。
- pooled DMR 来源：MethSCAn pooled 作业 `<job-id>`（run
  `Results/Methscan/<run-name>`）完成的 91/91 组比较。
- DMR 过滤：raw `p < 0.01` 与 `abs(methdiff) >= 0.25`。
- 不做 Top200 或其他 Top-N 截断。
- 精确去重 key：hypo cell type + chromosome + start + end。
- blacklist：移除重叠至少 20% 的源区间。
- DMR 合并：合并所有重叠的保留区间。
- VMR/DMR 重叠处理：保留已选 VMR，只追加零重叠的 DMR。

Preparation completed with 1,642,437 qualifying rows, 1,104,139 exact-unique
hypo-DMRs, and 321,910 merged intervals spanning 967,355,000 bp. The longest
merged interval is 58 kb.

准备步骤完成时共有 1,642,437 行符合条件、1,104,139 个完全唯一的 hypo-DMR、以及
321,910 个合并区间，总跨度 967,355,000 bp。最长的合并区间为 58 kb。

## Models and outputs / 模型与输出

Six models are defined as MethSCAn VMR `{0.01,0.02,0.05} x {Top10k,Top30k}`,
with the model-specific non-overlapping pooled DMR set appended. Each model uses
integer CGN `mc/cov`, batch key `sample_id`, seed 0, up to 500 epochs, and
validation-ELBO early stopping.

6 个模型定义为 MethSCAn VMR `{0.01,0.02,0.05} x {Top10k,Top30k}`，并追加各自
模型专属的非重叠 pooled DMR 集合。每个模型都使用整数 CGN `mc/cov`、batch key
`sample_id`、seed 0、最多 500 epochs，以及验证集 ELBO early stopping。

Each completed route must contain:

- `input.h5mu` and `input.summary.json`;
- saved MethylVI model and `methylvi_embedding.h5ad`;
- cell-type, sample, condition, and Leiden UMAP PNGs;
- supervised UMAP and methylation QC outputs;
- `input.COMPLETE` and `model.COMPLETE`.

每个已完成的路线必须包含：

- `input.h5mu` 与 `input.summary.json`；
- 保存的 MethylVI 模型与 `methylvi_embedding.h5ad`；
- cell-type、sample、condition 与 Leiden 的 UMAP PNG；
- supervised UMAP 与甲基化 QC 输出；
- `input.COMPLETE` 与 `model.COMPLETE`。

The supervised UMAP is exported at weights 0.2/0.5/0.7/0.9. Weights 0.5/0.7/0.9
distort latent structure and are sensitivity analysis only; `target_weight=0.2`
is the primary supervised view.

supervised UMAP 以权重 0.2/0.5/0.7/0.9 导出。权重 0.5/0.7/0.9 会扭曲 latent
结构，仅用于敏感性分析；`target_weight=0.2` 才是主要的 supervised 视图。

## Current status / 当前状态

Status date / 状态日期: `<date>`. Fill in the run's own jobs, resources, and
status; the workflow runs serially on a single node:

状态日期：`<date>`。按下表填写本次 run 的作业、资源与状态；流程在单节点上串行
执行：

| Stage / 阶段 | Job(s) / 作业 | Resources / 资源 | Status / 状态 |
|---|---|---|---|
| Pooled-DMR counts / pooled-DMR 计数 | `<job-id>` | 40 CPU / 48G | `<status>` |
| Six input joins / 6 个输入 join | `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>` | 4 CPU / 40G each / 各 4 CPU / 40G | `<status>` |
| Six MethylVI trains/plots / 6 个 MethylVI 训练与绘图 | `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>`, `<job-id>` | 32 CPU / 48G each / 各 32 CPU / 48G | `<status>` |
| Final validation / 最终校验 | `<job-id>` | 2 CPU / 8G | `<status>` |

Result root:

结果根目录：

`Results/MethylVI/<run-name>`.

Status is **in progress** until all six models pass. Only the final validation
job may write the final `run_summary.json` and the `workflow.COMPLETE` marker
after all six models pass.

在全部 6 个模型通过之前，状态为 **in progress**。只有在全部 6 个模型通过之后，
最终校验作业才可以写出最终的 `run_summary.json` 与 `workflow.COMPLETE` 标记。
