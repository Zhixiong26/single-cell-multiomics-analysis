# MethylVI workflows / MethylVI 流程

This directory contains three coordinated MethylVI feature routes over the same
`$SCMO_EXPECTED_CELLS` MethSCAn-filtered cells (reference project example:
6,264). The completed baseline experiment builds four feature inputs and eight
ALLCools/VMR models; the VMR+pooled-DMR extension adds six models.

本目录包含三条协同的 MethylVI 特征路线，均作用于同一批 `$SCMO_EXPECTED_CELLS`
个 MethSCAn filter 后的细胞（参考项目示例：6,264）。已完成基线实验构建 4 个
特征输入、8 个 ALLCools/VMR 模型；VMR+pooled-DMR 扩展再增加 6 个模型。

## Route layout / 路线布局

```text
Bismark coverage / ALLC
├── allcools/       ALLCools 5-kb feature selection → MethylVI
├── MethSCAn output VMRs
│   └── vmr/        selected MethSCAn VMR BED → MethylVI
└── pooled cell-type DMRs + selected VMRs
    └── vmr_dmr/    all unique pooled hypo-DMRs + VMR → MethylVI

shared/             training, supervised-UMAP, and methylation-QC code
```

`allcools/`、`vmr/`、`vmr_dmr/` 是三条特征路线；`shared/` 提供训练、
supervised UMAP 与甲基化 QC 代码。

## Production submission / 正式提交

```bash
bash Scripts/Methylvi/submit_methylvi_8models.sh
```

The dependency graph builds the shared feature inputs first, trains eight
models in parallel when their inputs are ready, and runs a final output
summary. The experiment root must not already exist.

依赖图先构建共享特征输入，输入就绪后并行训练 8 个模型，最后运行输出汇总。
实验根目录必须不存在（避免覆盖既有结果）。

The completed baseline is the 8-model run root `Results/MethylVI/<run-name>`
(8/8 models complete in the reference project). Read completion from the run's
manifest (`run_summary.tsv` / `run_summary.json` plus `workflow.COMPLETE`), not
from a dated directory name. The VMR+DMR extension is submitted separately:

基线结果是 8 模型 run 根目录 `Results/MethylVI/<run-name>`（参考项目中 8/8 模型
完成）。完成状态以 run manifest（`run_summary.tsv` / `run_summary.json` 与
`workflow.COMPLETE`）为准，而不是按带日期的目录名判断。VMR+DMR 扩展单独提交：

```bash
bash Scripts/Methylvi/vmr_dmr/submit_vmr_dmr_pipeline.sh
```

It writes only to its own `Results/MethylVI/<run-name>` root, reuses the
completed baseline VMR inputs, and serializes its six join/train pairs on a
single node.

它只写入自己的 `Results/MethylVI/<run-name>` 根目录，复用已完成的基线 VMR 输入，
并把 6 组 join/train 串行化在同一个节点上执行。

Individual route runners remain available for verification and diagnosis:

各路线 runner 保留用于校验与诊断：

```bash
bash Scripts/Methylvi/allcools/run.sh <stage>
bash Scripts/Methylvi/vmr/run.sh <stage>
```

Production execution is available only through the project DAG
(`tools/plan_workflow.py` and `tools/submit_workflow.py`); stage-local Slurm
submission is intentionally unsupported.

正式执行只能通过项目 DAG（`tools/plan_workflow.py` 与 `tools/submit_workflow.py`）；
不支持绕过证据链的阶段级 Slurm 直提。

## Cell-selection authority / 细胞名单权威

All routes use `03_filtered/column_header.txt` as the sole cell-selection
authority (the MethSCAn-filtered cell list). The selected-ALLC manifest maps
those cell IDs to original indexed ALLCs and is a path-resolution table, not a
second selection rule. The VMR route also requires a completed MethSCAn
`run_summary.json` before preparing counts. The upstream MethSCAn run is
selected by run name (`Results/Methscan/<run-name>`), not by a dated path.

所有路线以 `03_filtered/column_header.txt`（MethSCAn filter 后的细胞名单）作为
唯一细胞选择权威。selected-ALLC manifest 只把这些 cell ID 映射到原始 indexed
ALLC，属于路径解析表，不是第二条选择规则。VMR 路线在准备计数前还要求 MethSCAn
`run_summary.json` 已完成。上游 MethSCAn run 按 run 名选择
（`Results/Methscan/<run-name>`），不依赖带日期的路径。

Note / 注意：cell-type 标签属于逐路线配置，不是硬编码默认值。请显式设置
`SCMO_CELLTYPE_KEY`（默认 `cell_type`）；参考项目的 ALLCools 路线使用的是
`manual_celltype`，因此该值必须逐路线确认。含 `/` 的 cell-type 标签（例如
`Secretory / mucous epithelial`）不能作为 AnnData `uns` 的字典 key，需按
`shared/README.md` 的方式以两个对齐数组加 JSON sidecar 保存。

Note / 注意: set `SCMO_CELLTYPE_KEY` explicitly per route (default `cell_type`);
the reference project's ALLCools route used `manual_celltype`. Cell-type labels
containing `/` cannot be AnnData `uns` dictionary keys — store them as two
aligned arrays plus a JSON sidecar, as described in `shared/README.md`.

## Feature builders and the eight models / 特征构建与八个模型

Four feature builders feed eight models: ALLCools 5-kb × {Top10k,Top30k} and
MethSCAn VMR {0.01,0.02,0.05} × {Top10k,Top30k}. MethylVI uses integer CGN
`mc/cov`, trains with `sample_id` as batch key, and exports ordinary
UMAP/Leiden, supervised UMAP, sequencing-depth, overall-mCG, and mean-mCG
diagnostics.

4 个特征构建任务支撑 8 个模型：ALLCools 5-kb × {Top10k,Top30k}，以及 MethSCAn
VMR {0.01,0.02,0.05} × {Top10k,Top30k}。MethylVI 使用整数 CGN `mc/cov`，以
`sample_id` 作为 batch key 训练，并导出普通 UMAP/Leiden、supervised UMAP、
测序深度、整体 mCG 与平均 mCG 诊断。

## VMR+DMR definition / VMR+DMR 定义

The VMR+DMR route uses all pooled hypo-DMRs passing raw `p < 0.01` and
`abs(methdiff) >= 0.25`, without Top200/Top-N truncation. Exact intervals are
deduplicated within hypo cell type, overlapping DMRs are merged, and only DMRs
with zero overlap to a model's selected VMRs are appended. Its preparation
produced 1,104,139 exact-unique and 321,910 merged DMR intervals in the
reference project.

VMR+DMR 路线使用全部满足 raw `p < 0.01` 与 `abs(methdiff) >= 0.25` 的 pooled
hypo-DMR，不做 Top200/Top-N 截断。完全相同的区间在同一 hypo cell type 内去重，
重叠的 DMR 合并，且只追加与某个模型所选 VMR 零重叠的 DMR。参考项目中该准备
步骤产出 1,104,139 个完全唯一区间与 321,910 个合并区间。

The high supervised UMAP weights (0.5/0.7/0.9) distort latent structure and are
sensitivity analysis only; `target_weight=0.2` is the primary supervised view.

高 supervised UMAP 权重（0.5/0.7/0.9）会扭曲 latent 结构，仅用于敏感性分析；
`target_weight=0.2` 才是主要的 supervised 视图。

Read the corresponding workflow README before running. Override the run,
variance branch, or BED only with a completed MethSCAn result.

运行前请阅读对应的流程 README。只有在 MethSCAn 结果已完成时，才覆盖 run、
variance 分支或 BED。
