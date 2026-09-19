# ALLCools → MethylVI

This workflow starts from the `$SCMO_EXPECTED_CELLS` cells retained in the
completed MethSCAn filter header, resolves their original indexed ALLCs through
the MethSCAn input manifest, and uses ALLCools 5-kb hypo-score regions for
feature selection and clustering, then reconstructs integer CGN `mc/cov` counts
for MethylVI. The ALLCools score itself is never used as a MethylVI count
matrix.

本流程从已完成 MethSCAn filter header 中保留的 `$SCMO_EXPECTED_CELLS` 个细胞
出发，通过 MethSCAn input manifest 解析其原始 indexed ALLC，使用 ALLCools 5-kb
hypo-score 区域做特征选择与聚类，再为 MethylVI 重建整数 CGN `mc/cov` 计数。
ALLCools 的 score 本身永不作为 MethylVI 计数矩阵。

## Stages / 阶段

| Command / 命令 | Result / 结果 |
|---|---|
| `verify` | Validate input paths, environments, and cell inventory. / 校验输入路径、环境与细胞清单。 |
| `prepare` | Reuse or create ALLC files and generate the 5-kb MCDS. / 复用或创建 ALLC 文件并生成 5-kb MCDS。 |
| `cluster` | Rank blacklist-filtered 5-kb bins by current-cell hypo prevalence and select exact top 30k. / 将 blacklist 过滤后的 5-kb bins 按当前细胞 hypo 阳性率排序，并精确选出 top 30k。 |
| `build` | Aggregate retained bins to integer MethylVI `mc/cov` layers. / 将保留的 bins 聚合为整数 MethylVI `mc/cov` 层。 |
| `train` | Train MethylVI and save the latent embedding. / 训练 MethylVI 并保存 latent embedding。 |
| `plots-before` / `plots-after` / `supervised` | Export the respective UMAPs. / 分别导出对应的 UMAP。 |

Run from the repository root:

在项目根目录下运行：

```bash
bash Scripts/Methylvi/allcools/run.sh verify
bash Scripts/Methylvi/allcools/run.sh all
```

The production DAG builds 30k once and derives a strictly nested top-10k H5MU
using `selection_rank`; it does not rescan ALLCs. The effective
`MVI_HYPO_PERCENT` and hard feature-count check are written to
`feature_filter_summary.json`. Configure paths, feature selection,
environments, and training parameters in `00_methylvi_config.sh` or by exporting
the documented `SCMO_*` variables before invoking `run.sh`.

正式 DAG 只构建一次 30k，并用 `selection_rank` 派生严格嵌套的 top-10k H5MU，
不重复扫描 ALLC。生效的 `MVI_HYPO_PERCENT` 与硬性特征数校验会写入
`feature_filter_summary.json`。路径、特征选择、环境与训练参数可在
`00_methylvi_config.sh` 中配置，或在调用 `run.sh` 前导出文档化的 `SCMO_*`
变量。

Note / 注意：请逐路线显式设置 `SCMO_CELLTYPE_KEY`（默认 `cell_type`）。参考
项目中 ALLCools 路线使用的是 `manual_celltype`；该 key 决定聚类标签与
supervised UMAP 的目标列，配置错误会静默改变细胞类型视图。

Production execution is available only through the project DAG; the stage-local
Slurm wrappers under `slurm/` are not a supported production entry point.

正式执行只能通过项目 DAG；`slurm/` 下的阶段级 Slurm 封装不是受支持的正式入口。

The baseline experiment root (`Results/MethylVI/<run-name>`, 8/8 models) has
completed both ALLCools models (Top10k and Top30k), including embeddings,
ordinary and supervised UMAPs, methylation QC, and `model.COMPLETE` markers.

基线实验根目录（`Results/MethylVI/<run-name>`，8/8 模型）中两个 ALLCools 模型
（Top10k 与 Top30k）均已完成，包括 embedding、普通与 supervised UMAP、甲基化
QC 以及 `model.COMPLETE` 标记。
