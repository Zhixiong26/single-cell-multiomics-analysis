# MethylVI report / MethylVI 报告

## Maintained workflow contract / 当前流程约定

All production routes use the same `$SCMO_EXPECTED_CELLS` MethSCAn-filtered
cells (reference project: 6,264 = CYL 2,919 + ZCP 3,345; sample names are
examples only), original indexed ALLCs, integer CGN methylated/coverage
(`mc/cov`) counts, and `sample_id` as the MethylVI batch key. Transformed
methylation scores and heatmap ratios are never used as MethylVI count
matrices.

所有正式路线统一使用 MethSCAn filter 后的 `$SCMO_EXPECTED_CELLS` 个细胞（参考
项目为 6,264，即 CYL 2,919 + ZCP 3,345；样本名仅为示例）、原始 indexed ALLC、
整数 CGN `mc/cov` 计数，并以 `sample_id` 作为 batch key。任何变换后的 score 或
热图 ratio 都不会作为 MethylVI 训练计数。

The repository maintains three routes:

1. `allcools/`: ALLCools-selected 5-kb bins, Top10k/Top30k.
2. `vmr/`: MethSCAn VMR thresholds 0.01/0.02/0.05, each Top10k/Top30k.
3. `vmr_dmr/`: each selected VMR input plus all qualifying, merged,
   non-overlapping pooled cell-type DMRs.

本仓库维护三条路线：

1. `allcools/`：ALLCools 选出的 5-kb bins，Top10k/Top30k。
2. `vmr/`：MethSCAn VMR 阈值 0.01/0.02/0.05，各自 Top10k/Top30k。
3. `vmr_dmr/`：每个已选 VMR 输入，加上全部符合条件、已合并、互不重叠的 pooled
   cell-type DMR。

## Completed 10k/30k baseline / 已完成基线

The baseline 8-model run root `Results/MethylVI/<run-name>` is complete
(`workflow.COMPLETE`; final summary job `<job-id>`). Its experiment-level
`run_summary.tsv` records 8/8 models as complete:

- ALLCools 5-kb: Top10k and Top30k;
- MethSCAn VMR: `{0.01,0.02,0.05} x {Top10k,Top30k}`.

基线 8 模型 run 根目录 `Results/MethylVI/<run-name>` 已完成（`workflow.COMPLETE`；
最终汇总作业 `<job-id>`）。其实验级 `run_summary.tsv` 记录 8/8 模型完成：

- ALLCools 5-kb：Top10k 与 Top30k；
- MethSCAn VMR：`{0.01,0.02,0.05} x {Top10k,Top30k}`。

Every model contains a saved MethylVI model, 20-dimensional latent embedding,
ordinary UMAPs colored by cell type/sample/condition/Leiden cluster,
supervised UMAPs, and methylation QC. The VMR source counts were 39,553,
80,818, and 166,618 before the MethylVI route's canonical/blacklist/coverage
selection.

每个模型都包含保存的 MethylVI 模型、20 维 latent embedding、按
cell type/sample/condition/Leiden cluster 着色的普通 UMAP、supervised UMAP 以及
甲基化 QC。在经过 MethylVI 路线的 canonical/blacklist/coverage 选择之前，VMR
源计数分别为 39,553、80,818、166,618。

## Train/validation split / 训练验证划分

The shared trainer does not define an independent test set. With scvi-tools
1.5.0.post1 defaults and seed 0, all cells are globally shuffled (not
stratified by sample or cell type). For `$SCMO_EXPECTED_CELLS` = 6,264
(reference project):

- training: 5,638 cells (CYL 2,638; ZCP 3,000);
- validation: 626 cells (CYL 281; ZCP 345);
- test: 0 cells.

共享训练器不定义独立测试集。在 scvi-tools 1.5.0.post1 默认设置与 seed 0 下，
所有细胞全局打乱（不按 sample 或 cell type 分层）。当 `$SCMO_EXPECTED_CELLS`
= 6,264（参考项目）时：

- 训练：5,638 个细胞（CYL 2,638；ZCP 3,000）；
- 验证：626 个细胞（CYL 281；ZCP 345）；
- 测试：0 个细胞。

Validation ELBO controls early stopping. Final latent embeddings and UMAPs are
computed for all cells, so the UMAP is not an independent test-set performance
estimate.

验证集 ELBO 控制 early stopping。最终 latent embedding 与 UMAP 对全部细胞计算，
因此 UMAP 不是独立的测试集性能估计。

## Pooled DMR upstream / pooled DMR 上游

MethSCAn pooled job `<job-id>` completed on `<date>`. It combined cells from
both samples with the same RNA cell type and completed all 91 pairwise
comparisons among 14 cell types using all `$SCMO_EXPECTED_CELLS` cells. There
were 0 failures and 0 fallback comparisons, with 6,437,767 output DMR rows.

MethSCAn pooled 作业 `<job-id>` 于 `<date>` 完成。它按相同 RNA cell type 合并
两个样本的细胞，使用全部 `$SCMO_EXPECTED_CELLS` 个细胞完成 14 个 cell type 之间
的全部 91 组两两比较。0 个失败、0 个 fallback 比较，输出 6,437,767 行 DMR。

For the VMR+DMR route, every pooled DMR row satisfying raw `p < 0.01` and
`abs(methdiff) >= 0.25` is retained; there is no Top200 or other Top-N
truncation. After the shared 20% blacklist rule:

- 1,642,437 rows passed the thresholds before exact deduplication;
- 1,104,139 DMRs were unique by hypo cell type and exact coordinates;
- overlap merging produced 321,910 non-overlapping DMR intervals;
- merged span: 967,355,000 bp; maximum interval length: 58 kb.

对 VMR+DMR 路线，所有满足 raw `p < 0.01` 与 `abs(methdiff) >= 0.25` 的 pooled
DMR 行全部保留，不做 Top200 或其他 Top-N 截断。经过共享的 20% blacklist 规则后：

- 精确去重前有 1,642,437 行通过阈值；
- 1,104,139 个 DMR 按 hypo cell type 与精确坐标唯一；
- 重叠合并后得到 321,910 个互不重叠的 DMR 区间；
- 合并跨度 967,355,000 bp；最大区间长度 58 kb。

For each model, the selected Top10k/Top30k VMR intervals remain unchanged and
only pooled DMR intervals with zero genomic overlap to those VMRs are appended,
avoiding duplicate representation of the same CpGs.

对每个模型，已选的 Top10k/Top30k VMR 区间保持不变，只追加与这些 VMR 零基因组
重叠的 pooled DMR 区间，避免同一批 CpG 被重复表示。

## VMR+DMR run status / VMR+DMR 运行状态

Status date / 状态日期: `<date>`. Results are isolated under the VMR+DMR run root
`Results/MethylVI/<run-name>`; the completed baseline is not overwritten. Record
the run's own jobs, resources, and status here:

结果隔离在 VMR+DMR run 根目录 `Results/MethylVI/<run-name>` 下；已完成的基线不会
被覆盖。在此记录本次 run 的作业、资源与状态：

| Stage / 阶段 | Job(s) / 作业 | Resources / 资源 | Status / 状态 |
|---|---|---|---|
| Build pooled-DMR integer counts / 构建 pooled-DMR 整数计数 | `<job-id>` | 40 CPU / 48G | `<status>` |
| var0.01 Top10k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.01 Top10k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| var0.01 Top30k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.01 Top30k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| var0.02 Top10k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.02 Top10k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| var0.02 Top30k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.02 Top30k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| var0.05 Top10k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.05 Top10k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| var0.05 Top30k join / join | `<job-id>` | 4 CPU / 40G | `<status>` |
| var0.05 Top30k train / 训练 | `<job-id>` | 32 CPU / 48G | `<status>` |
| Final validation and summary / 最终校验与汇总 | `<job-id>` | 2 CPU / 8G | `<status>` |

The six join/train pairs are serialized. Join jobs request 4 CPU/40G; training
jobs request 32 CPU/48G; summary requests 2 CPU/8G. Each trainer automatically
exports standard UMAPs for `cell_type`, `sample_id`, `condition`, and
`methylVI_leiden`, followed by supervised UMAP and methylation QC. Only the
final job may create `workflow.COMPLETE`; pending jobs must not be described as
complete before that marker and `run_summary.json` exist.

6 组 join/train 串行执行。join 作业申请 4 CPU/40G；训练作业申请 32 CPU/48G；
汇总作业申请 2 CPU/8G。每个 trainer 自动导出 `cell_type`、`sample_id`、
`condition` 与 `methylVI_leiden` 的标准 UMAP，随后是 supervised UMAP 与甲基化
QC。只有最终作业可以创建 `workflow.COMPLETE`；在该标记与 `run_summary.json`
存在之前，等待中的作业不得描述为已完成。
