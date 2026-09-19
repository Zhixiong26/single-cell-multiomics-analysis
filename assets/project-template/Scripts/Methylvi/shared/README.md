# Shared MethylVI steps / 共享 MethylVI 步骤

`04_train_methylvi.py`, `05_plot_supervised_umap.py`, and
`07_plot_methylation_qc.py` are shared by the ALLCools, MethSCAn-VMR, and
VMR+pooled-DMR workflows; `06_subset_methylvi_input.py` derives the nested
top-N subsets from `selection_rank`. They consume common `SCMO_*` MethylVI
environment variables, which each workflow runner sets for its own result
route.

`04_train_methylvi.py`、`05_plot_supervised_umap.py` 与
`07_plot_methylation_qc.py` 被 ALLCools、MethSCAn-VMR 和 VMR+pooled-DMR 三条
流程共用；`06_subset_methylvi_input.py` 依据 `selection_rank` 派生嵌套的 top-N
子集。它们读取统一的 `SCMO_*` MethylVI 环境变量，由各流程 runner 为自己的
结果路线设置。

Training uses the scvi-tools default random cell split with seed 0: for
`$SCMO_EXPECTED_CELLS` = 6,264 (reference project) this is 5,638 training cells,
626 validation cells, and no independent test set. The validation ELBO drives
early stopping; latent embeddings and UMAPs are subsequently generated for all
6,264 cells. Training uses up to 500 epochs and `batch_size=32`, with
`sample_id` as the batch key.

训练使用 scvi-tools 默认随机细胞划分与 seed 0：当 `$SCMO_EXPECTED_CELLS` =
6,264（参考项目）时为 5,638 个训练细胞、626 个验证细胞，没有独立测试集。验证集
ELBO 驱动 early stopping；随后对全部 6,264 个细胞生成 latent embedding 与
UMAP。训练最多 500 epochs、`batch_size=32`，并以 `sample_id` 作为 batch key。

Supervised UMAP is exported at `target_weight` 0.2/0.5/0.7/0.9. The high weights
distort latent structure and are sensitivity analysis only; `target_weight=0.2`
is the primary supervised view.

Supervised UMAP 以 `target_weight` 0.2/0.5/0.7/0.9 导出。高权重会扭曲 latent
结构，仅用于敏感性分析；`target_weight=0.2` 才是主要的 supervised 视图。

## Cell-type labels containing `/` / 含 `/` 的细胞类型标签

Cell-type labels may legitimately contain `/` (for example
`Secretory / mucous epithelial`). AnnData serializes `uns` mapping keys as HDF5
group names, so such a label cannot be used as a dictionary key in `uns` — the
write fails or the key is mangled. The supervised-UMAP step therefore stores the
label/code mapping as two aligned arrays (`target_labels` and `target_codes`)
inside `uns`, and writes the human-readable dictionary separately as JSON
(`target_mapping` in `supervised_umap_summary.json`). Do not reintroduce a
label-keyed `uns` dict; any reader must rebuild the mapping by zipping the two
aligned arrays.

细胞类型标签中合法地可能含有 `/`（例如 `Secretory / mucous epithelial`）。
AnnData 会把 `uns` 的 mapping key 序列化为 HDF5 group 名，因此这类标签不能直接
作为 `uns` 的字典 key——写入会失败或 key 被改写。因此 supervised UMAP 步骤把
label/code 映射保存为 `uns` 中两个对齐数组（`target_labels` 与
`target_codes`），并把人类可读字典单独写成 JSON（`supervised_umap_summary.json`
中的 `target_mapping`）。不要重新引入以 label 为 key 的 `uns` 字典；读取方必须
用两个对齐数组重建映射。
