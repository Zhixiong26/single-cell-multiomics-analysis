# MethSCAn VMRs → MethylVI / MethSCAn VMR → MethylVI

This workflow trains MethylVI on integer CGN `mc/cov` counts aggregated over a
VMR BED produced by a completed MethSCAn run. It is independent from the
ALLCools 5-kb workflow and writes to a separate result root.

本流程在已完成 MethSCAn run 产出的 VMR BED 上聚合整数 CGN `mc/cov` 计数来训练
MethylVI。它与 ALLCools 5-kb 流程相互独立，并写入独立的结果根目录。

## Choose the MethSCAn VMR set / 选择 MethSCAn VMR 集合

MethSCAn produces one VMR BED per variance threshold. Defaults target the
current formal run and its `0.01` branch (`Results/Methscan/<run-name>/04_scan/
var_<threshold>/VMRs.bed`). To use another completed run or branch, override the
run, the threshold, or the BED directly; for example:

MethSCAn 每个 variance 阈值产出一个 VMR BED。默认指向当前正式 run 及其 `0.01`
分支（`Results/Methscan/<run-name>/04_scan/var_<threshold>/VMRs.bed`）。如需改用
其他已完成的 run 或分支，可覆盖 run、阈值，或直接覆盖 BED，例如：

```bash
export SCMO_VMR_SOURCE_BED="$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>/04_scan/var_0.01/VMRs.bed"
```

The packaged route config exports this override as `VMR_SOURCE_BED`; set the
neutral `SCMO_VMR_SOURCE_BED` name where the project configuration maps it.

打包的路线配置以 `VMR_SOURCE_BED` 之名导出该覆盖项；在项目配置提供映射的位置请
使用中性的 `SCMO_VMR_SOURCE_BED` 名称。

`run.sh verify` refuses to proceed when the selected BED is missing. The VMR
builder uses `03_filtered/column_header.txt` as the sole cell list and uses the
selected-ALLC manifest only to resolve paths, so MethylVI uses the same
`$SCMO_EXPECTED_CELLS` post-filter cells (reference project: 6,264). VMRs are
restricted to canonical chromosomes, blacklist-filtered, checked for overlap,
and assigned the highest covered-cell cutoff that still permits the 30k target.
The resulting cell cutoff and effective percentage are recorded in
`build_summary.json`.

当选定的 BED 缺失时，`run.sh verify` 会拒绝继续。VMR 构建器以
`03_filtered/column_header.txt` 作为唯一细胞名单，selected-ALLC manifest 只用于
解析路径，因此 MethylVI 使用同一批 `$SCMO_EXPECTED_CELLS` 个 filter 后细胞
（参考项目：6,264）。VMR 被限制在 canonical 染色体上，经 blacklist 过滤、重叠
检查，并赋予仍能满足 30k 目标的最高 covered-cell cutoff。最终的 cell cutoff 与
生效百分比记录在 `build_summary.json`。

After coverage eligibility, VMRs are ranked by `peak_var` descending, then
`n_obs_cells`, `n_cpg`, and stable VMR ID. Each threshold builds exact top 30k
once and derives a nested top 10k without rescanning ALLCs.

通过 coverage 资格筛选后，VMR 依次按 `peak_var` 降序、`n_obs_cells`、`n_cpg`
和稳定 VMR ID 排序。每个阈值只构建一次精确 top 30k，并派生嵌套的 top 10k，
不重复扫描 ALLC。

Blacklist handling uses the shared rule (reference project: ENCODE GRCh38
blacklist, `blacklist_fraction=0.2`, MD5
`393688b4f06c9ce26165d47433dd8c37`, canonical chromosomes only).

blacklist 处理沿用共享规则（参考项目：ENCODE GRCh38 blacklist，
`blacklist_fraction=0.2`，MD5 `393688b4f06c9ce26165d47433dd8c37`，仅 canonical
染色体）。

## Stages / 阶段

```bash
bash Scripts/Methylvi/vmr/run.sh verify
bash Scripts/Methylvi/vmr/run.sh prepare
bash Scripts/Methylvi/vmr/run.sh build
bash Scripts/Methylvi/vmr/run.sh train
bash Scripts/Methylvi/vmr/run.sh plots
bash Scripts/Methylvi/vmr/run.sh all
```

Run all three thresholds at both feature targets through the project-level DAG:

通过项目级 DAG 运行三个阈值 × 两个特征目标：

```bash
bash Scripts/Methylvi/submit_methylvi_8models.sh
```

The route's output root is `Results/MethylVI/<run-name>`; the VMR route root
must be separate from the ALLCools route root. Paths, the variance branch, and
MethylVI parameters are in `00_vmr_methylvi_config.sh`. Production submission
goes through the project DAG; the stage-local Slurm wrappers under `slurm/` are
not a supported production entry point.

该路线输出根目录为 `Results/MethylVI/<run-name>`；VMR 路线根目录必须与 ALLCools
路线根目录分开。路径、variance 分支与 MethylVI 参数位于
`00_vmr_methylvi_config.sh`。正式提交通过项目 DAG 进行；`slurm/` 下的阶段级
Slurm 封装不是受支持的正式入口。

## Formal completed experiment / 已完成的正式实验

The baseline experiment root (`Results/MethylVI/<run-name>`, 8/8 models)
contains six completed VMR models: `{0.01,0.02,0.05} x {Top10k,Top30k}`. Every
route has `model.COMPLETE`, a 20-dimensional latent embedding, and UMAPs colored
by cell type, sample, condition, and MethylVI Leiden cluster. The same selected
VMR inputs are reused by `../vmr_dmr/`, which appends only non-overlapping pooled
DMR features and writes to a separate result root.

基线实验根目录（`Results/MethylVI/<run-name>`，8/8 模型）包含 6 个已完成的 VMR
模型：`{0.01,0.02,0.05} x {Top10k,Top30k}`。每条路线都有 `model.COMPLETE`、20 维
latent embedding，以及按 cell type、sample、condition 与 MethylVI Leiden
cluster 着色的 UMAP。同一批已选 VMR 输入会被 `../vmr_dmr/` 复用；后者只追加
互不重叠的 pooled DMR 特征，并写入独立的结果根目录。

Supervised UMAP weights 0.5/0.7/0.9 are sensitivity analysis only and distort
latent structure; `target_weight=0.2` is the primary supervised view. Training
uses integer `mc/cov`, `sample_id` batch key, seed 0, up to 500 epochs,
`batch_size=32`, and validation-ELBO early stopping.

supervised UMAP 权重 0.5/0.7/0.9 仅为敏感性分析，会扭曲 latent 结构；
`target_weight=0.2` 才是主要的 supervised 视图。训练使用整数 `mc/cov`、
`sample_id` batch key、seed 0、最多 500 epochs、`batch_size=32`，并由验证集
ELBO 控制 early stopping。
