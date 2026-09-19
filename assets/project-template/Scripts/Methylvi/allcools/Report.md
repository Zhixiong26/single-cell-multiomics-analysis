# ALLCools 5-kb → MethylVI report / ALLCools 5-kb → MethylVI 报告

## Maintained contract / 当前约定

The maintained cell list is the `$SCMO_EXPECTED_CELLS` IDs in the completed
MethSCAn `03_filtered/column_header.txt` (reference project: 6,264 = CYL 2,919 +
ZCP 3,345; sample names are examples only). The upstream selected manifest is
used only to resolve each ID to its original indexed ALLC. The historical
6,554-cell legacy route is fallback code, not the production cell-selection
rule.

当前唯一细胞名单是 MethSCAn `03_filtered/column_header.txt` 中的
`$SCMO_EXPECTED_CELLS` 个细胞（参考项目：6,264，即 CYL 2,919 + ZCP 3,345；样本名
仅为示例）。上游 manifest 只负责将 cell ID 映射到原始 indexed ALLC；历史
6,554-cell 遗留路线仅为兼容 fallback，不是正式细胞选择规则。

After 5-kb CGN MCDS generation and blacklist filtering, bins are ranked by the
number of current cells passing the binarized hypo-score cutoff. The production
branch selects exactly top 30,000 bins, records the effective
`MVI_HYPO_PERCENT` and boundary ties in `feature_filter_summary.json`, and
derives a strictly nested top-10,000 H5MU using `selection_rank`.

5-kb CGN MCDS 生成与 blacklist 过滤后，按当前细胞集中 binarized hypo-score 阳性
细胞数对 bins 排序。正式流程精确选择 top 30,000，并把生效的
`MVI_HYPO_PERCENT` 与边界并列情况记录到 `feature_filter_summary.json`；top
10,000 使用 `selection_rank` 派生，严格嵌套且不重复扫描 ALLC。

The reference project applied the ENCODE GRCh38 blacklist with
`blacklist_fraction=0.2` (MD5 `393688b4f06c9ce26165d47433dd8c37`), 5-kb bins,
and `mc_context=CGN`; these values come from configuration, not from this
document.

参考项目使用 ENCODE GRCh38 blacklist，`blacklist_fraction=0.2`
（MD5 `393688b4f06c9ce26165d47433dd8c37`）、5-kb bins 与 `mc_context=CGN`；这些
取值来自配置，而不是本文档。

The production experiment (`Results/MethylVI/<run-name>`, the 8-model baseline
run root) completed both ALLCools models (10k and 30k). Both have validated
integer-count inputs, `model.COMPLETE`, 20-dimensional embeddings, ordinary
cell-type/sample/condition/Leiden UMAPs, supervised UMAPs, and methylation QC
outputs. They are listed as `complete` in the experiment-level
`run_summary.tsv`.

正式实验（`Results/MethylVI/<run-name>`，8 模型基线 run 根目录）中的 ALLCools
Top10k 与 Top30k 两个模型均已完成，并通过输入、模型标记、embedding、UMAP 与
甲基化 QC 校验；实验级 `run_summary.tsv` 将两者均记录为 `complete`。

Note / 注意：ALLCools 路线的细胞类型列需要逐路线显式配置 `SCMO_CELLTYPE_KEY`
（默认 `cell_type`）；参考项目的 ALLCools 路线实际使用 `manual_celltype`。

Supervised UMAP views at weights 0.5/0.7/0.9 are sensitivity analyses only and
distort latent structure; `target_weight=0.2` is the primary supervised view.

权重 0.5/0.7/0.9 的 supervised UMAP 仅为敏感性分析，会扭曲 latent 结构；
`target_weight=0.2` 才是主要的 supervised 视图。
