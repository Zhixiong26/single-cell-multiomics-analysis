# MethSCAn workflow / MethSCAn 流程（整理版）

This directory is numbered in execution order and resolves the filename conflict between the two `02` scripts of the previous version. / 本目录按执行顺序编号，解决旧版两个 `02` 文件名冲突：

```text
00_methscan_config.sh
01_select_scanpy_cells.py
02_convert_allc_to_cov.py
03_prepare_methscan.py
04_vmr_scanpy.py
05_summarize_run.py
06_run_with_resources.py
07_methdiff_celltype.py
08_hypo_dmr_heatmaps.py
run_methscan.sbatch
run_methscan_common.sbatch
run_methscan_qc_stage.sbatch
run_methscan_branch.sbatch
run_methscan_methdiff.sbatch
run_methscan_hypo_heatmaps.sbatch
run_methscan_pooled_methdiff.sbatch
run_methscan_summary.sbatch
submit_methscan_pipeline.sh
submit_methscan_pooled_dmr.sh
```

A single submission runs selection, ALLC→COV, prepare, filter, and smooth in sequence; after smooth, pairwise `methscan diff` is executed for all retained cell types within each sample (example samples: `CYL`, `ZCP`), followed by raw-p fallback, Top200 hypo-DMR, the single-cell matrix, and the three heatmap sets. Three VMR threshold branches run in parallel at the same time: `0.01`, `0.02`, `0.05`. ALLC→COV conversion uses `ProcessPoolExecutor`; the worker count defaults to `SLURM_CPUS_PER_TASK` and can be overridden with `SCMO_COV_WORKERS`. / 一次提交会依次运行 selection、ALLC→COV、prepare、filter、smooth；smooth 后在各样本内（示例/example: `CYL`、`ZCP`）执行所有保留 cell type 的两两 `methscan diff`，再完成 raw-p fallback、Top200 hypo-DMR、单细胞矩阵和三套热图。同时并行完成三个 VMR threshold 分支：`0.01`、`0.02`、`0.05`。ALLC→COV 转换使用 `ProcessPoolExecutor`，worker 默认自动等于 `SLURM_CPUS_PER_TASK`，也可用 `SCMO_COV_WORKERS` 覆盖。

The hypo-DMR downstream convention follows the workflow recorded by the reference project (example/示例: `scLC_ICI_PBMC`, commit `0c5c46c101fd669edd8b5c6e4e7dccb17d1f4bc3`). In that project the later-added "08 raw-p fallback" sits logically between pairwise diff and the original "05 Top200 / 06 matrix"; here the stages are written in true dependency order as `08_rawp_fallback`, `09_top200_hypo_DMRs`, `10_single_cell_DMR_matrix`, `11_hypo_DMR_heatmaps`. / hypo-DMR 下游口径借鉴参考项目（示例/example: `scLC_ICI_PBMC`，提交 `0c5c46c101fd669edd8b5c6e4e7dccb17d1f4bc3`）所记录的流程。该项目中后补的“08 raw-p fallback”在逻辑上位于 pairwise diff 与原“05 Top200 / 06 matrix”之间；本项目按真实依赖顺序统一写入 `08_rawp_fallback`、`09_top200_hypo_DMRs`、`10_single_cell_DMR_matrix`、`11_hypo_DMR_heatmaps`。

Use the one-click submitter to compose the common stages, the three threshold branches, and the final summary into a Slurm dependency DAG; the three branches run in parallel once the common stages complete. The original `run_methscan.sbatch` is still retained as a single-job compatible entry point. / 推荐使用一键提交器，将公共阶段、三个 threshold 分支和最终汇总组成 Slurm 依赖 DAG；三个分支在公共阶段完成后并行运行。原 `run_methscan.sbatch` 仍保留为单作业兼容入口。

## Run / 运行

ALLC must already be decompressed into the configured ALLC source directory (`$SCMO_ALLC_SOURCES`; example/示例: `Data/ALLCools`); the scripts never decompress, copy, or modify source files. The output directory must be a new directory that does not yet exist. / ALLC 必须已解压到配置的 ALLC 源目录（`$SCMO_ALLC_SOURCES`；示例/example: `Data/ALLCools`）；脚本不会解压、复制或修改源文件。输出目录必须是不存在的新目录：

```bash
sbatch <scheduler-profile-resource-options> \
  "$SCMO_PROJECT_ROOT/Scripts/Methscan/run_methscan.sbatch" \
  "$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>"
```

Partition, node pinning, CPU, and memory come from the scheduler profile, not from hard-coded values. The stage-by-stage resource shape below is an illustrative profile example only; the configured run root is `$SCMO_METHSCAN_RUN_DIR`. / 分区、节点、CPU 与内存由 scheduler profile 决定，不写死在命令里。下面的分阶段资源形状仅为示例 profile；配置的 run 根目录为 `$SCMO_METHSCAN_RUN_DIR`。

| Stage / 阶段 | Illustrative profile / 示例资源形状 |
|---|---|
| selection + ALLC→COV（公共作业 / common job） | 55 CPU / 250G，节点不固定 / node not pinned |
| prepare、filter、smooth（各一个串行作业 / one job each） | 4 CPU / 16G |
| 样本内 pairwise meth-diff / per-sample pairwise meth-diff | 8 CPU / 16G（2 comparisons × 4 threads） |
| hypo-DMR 矩阵与热图 / hypo-DMR matrix and heatmaps | 8 workers / 8 CPU / 48G |
| 每个 VMR threshold 分支 / each VMR threshold branch | 18 CPU / 80G |

Parallel branch submission / 并行分支提交：

```bash
"$SCMO_PROJECT_ROOT/Scripts/Methscan/submit_methscan_pipeline.sh" \
  "$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>"
```

This command submits the common job, prepare, filter, smooth, pairwise meth-diff, hypo-DMR/matrix/plotting, the 3 threshold branches, and the summary job; stages are chained automatically through `afterok`. / 该命令提交公共作业、prepare、filter、smooth、pairwise meth-diff、hypo-DMR/矩阵/绘图、3 个 threshold 分支和汇总作业；各阶段通过 `afterok` 自动串联。

For run directories whose common filtering and smoothing are already complete, the merged-cell-type DMR branch can be submitted separately / 对已经完成公共过滤和平滑的运行目录，可另外提交合并细胞类型 DMR 分支：

```bash
"$SCMO_PROJECT_ROOT/Scripts/Methscan/submit_methscan_pooled_dmr.sh" \
  "$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>"
```

This branch merges cells of the same-named RNA cell type across the samples into the pooled label (example/示例: `CYL` + `ZCP` → `CYL_ZCP`), then runs pairwise `methscan diff` over the cell types (example/示例: 14 cell types). It computes DMRs only; results are written independently to `<run>/pooled/07_methdiff/` and never overwrite the per-sample DMRs, and it does not run Top200, the single-cell matrix, or the heatmaps. / 该分支将各样本中同名 RNA cell type 的细胞合并为 pooled 标签（示例/example: `CYL`、`ZCP` 合并为 `CYL_ZCP`），然后进行细胞类型的两两 `methscan diff`（示例/example: 14 种细胞类型）。它只计算 DMR，结果独立写入 `<run>/pooled/07_methdiff/`，不会覆盖原有分样本 DMR，也不会运行 Top200、单细胞矩阵或热图。

The production pooled job `<job-id>` has completed: `<N>` filtered cells (example/示例: 6,264), `<N>` cell types (example/示例: 14), `<N>`/`<N>` successful comparisons (example/示例: 91/91), `<N>` failures and `<N>` fallbacks (example/示例: 0 and 0), and `<N>` DMR rows in total (example/示例: 6,437,767). This pooled DMR is later consumed by `Scripts/Methylvi/vmr_dmr/` under the all-unique-DMR convention of `raw p < 0.01` and `|methdiff| >= 0.25` to build the VMR+DMR MethylVI input; this does not change this directory's original per-sample Top200 heatmap convention. / 正式 pooled 作业 `<job-id>` 已完成：`<N>` 个过滤后细胞（示例/example: 6,264）、`<N>` 种细胞类型（示例/example: 14）、`<N>`/`<N>` 个比较成功（示例/example: 91/91）、`<N>` 个失败与 `<N>` 个 fallback（示例/example: 0 与 0），共输出 `<N>` 行 DMR（示例/example: 6,437,767）。该 pooled DMR 后续由 `Scripts/Methylvi/vmr_dmr/` 使用 `raw p < 0.01`、`|methdiff| >= 0.25` 的全量 unique-DMR 口径构建 VMR+DMR MethylVI 输入；这不会改变本目录原有的分样本 Top200 热图口径。

Default inputs are the configured ALLC source directory (example/示例: `Data/ALLCools`) and the configured RNA annotation table (`$SCMO_ANNOTATION`; example/示例: `Results/Scanpy/<notebook>/cell_id_cell_type.tsv`). Paths and parameters can be overridden through the environment variables in `00_methscan_config.sh`. / 默认输入为配置的 ALLC 源目录（示例/example: `Data/ALLCools`）与配置的 RNA 注释表（`$SCMO_ANNOTATION`；示例/example: `Results/Scanpy/<notebook>/cell_id_cell_type.tsv`）。可通过 `00_methscan_config.sh` 中的环境变量覆盖路径或参数。

## Execution stages / 执行阶段

1. `01_select_scanpy_cells.py`: discovers ALLC for the configured samples (example/示例: `CYL`, `ZCP`) and matches them against the RNA `cell_id`; `NA`, empty, and unannotated cells are excluded. / `01_select_scanpy_cells.py`：发现样本 ALLC（示例/example: `CYL`、`ZCP`），与 RNA `cell_id` 匹配；排除 `NA`、空值和未注释细胞。
2. `02_convert_allc_to_cov.py`: converts the selected ALLC into CpG-only Bismark `.cov.gz` and writes conversion QC. / `02_convert_allc_to_cov.py`：将入选 ALLC 转为 CpG-only Bismark `.cov.gz`，写入 conversion QC。
3. `03_prepare_methscan.py`: calls `methscan prepare --input-format bismark`. / `03_prepare_methscan.py`：调用 `methscan prepare --input-format bismark`。
4. `methscan filter`, `smooth`: use the default QC and smoothing parameters. / `methscan filter`、`smooth`：使用默认 QC 和平滑参数。
5. `07_methdiff_celltype.py`: first restricts to the chromosomes declared by `references.chrom_sizes` with a hard-link view that copies no data (the example project declares chr1–22/X/Y; the stage itself takes whatever that file lists and the summary records the scope it used); by default it runs pairwise `methscan diff` between RNA cell types within each sample, and merges same-named cell types across samples when `--pool-samples-as <pooled-label>` is passed (example/示例: `CYL_ZCP`). Defaults are `min-cells=6`, bandwidth=2,000, stepsize=1,000, threshold=0.02. The production wrapper enables `--resume`, which only reuses comparisons whose identity fields match, whose status is complete, and whose DMR format is valid. / `07_methdiff_celltype.py`：先以无数据复制的 hard-link view 限定 `references.chrom_sizes` 声明的染色体（示例项目声明 chr1–22/X/Y；该阶段本身按该文件列出的内容执行，并在 summary 里记录实际范围）；默认在每个样本内进行 RNA cell type 两两 `methscan diff`，传入 `--pool-samples-as <pooled-label>` 时则跨样本合并同名细胞类型（示例/example: `CYL_ZCP`）。默认 `min-cells=6`、bandwidth=2,000、stepsize=1,000、threshold=0.02。正式包装器启用 `--resume`，仅复用身份字段一致、状态完成且 DMR 格式有效的 comparison。
6. `08_hypo_dmr_heatmaps.py`: performs an isolated raw-p fallback only for the known `calc_fdr` division-by-zero failure; determines the hypo cell type from column 10, keeps DMRs with raw p `<0.01` and between-group mean ratio difference `>=0.25`, and takes the Top200 per cell type in descending difference order before merging overlapping intervals. / `08_hypo_dmr_heatmaps.py`：仅对已知 `calc_fdr` 除零错误执行隔离的 raw-p fallback；按第 10 列确定 hypo 细胞类型，保留 raw p `<0.01`、组间 mean ratio 差值 `>=0.25` 的 DMR，并按差值降序为每种细胞取 Top200 后合并重叠区间。
7. The same script computes the single-cell × merged-DMR matrix from `01_cov/cov/*.cov.gz`. Each value is the equal-weight arithmetic mean of the unique CpG ratio inside the DMR; missing coverage is `NA`. / 同一脚本从 `01_cov/cov/*.cov.gz` 计算单细胞×合并 DMR 矩阵。每个值是 DMR 内 unique CpG ratio 的等权算术平均，缺覆盖为 `NA`。
8. Per sample, only cell types assigned at least one of their own hypo-DMRs are kept as heatmap rows, producing `mean_ratio.png`, DMR-wise `zscore.png`, and `zscore_colorbar_compressed.png`. The compressed-colorbar figure keeps the original z-score values and only sets the color range to `[-1,1]`; the ordinary z-score figure displays the range `[-3,3]` by default. / 每个样本只保留至少分配到一个自身 hypo-DMR 的 cell type 作为热图行，生成 `mean_ratio.png`、DMR-wise `zscore.png` 和 `zscore_colorbar_compressed.png`。压缩 colorbar 图保留原 z-score 数值，只将颜色范围设为 `[-1,1]`；普通 z-score 图显示范围默认 `[-3,3]`。
9. `methscan scan → matrix → 04_vmr_scanpy.py`: run separately for each of the three VMR threshold branches, `0.01`/`0.02`/`0.05`. / `methscan scan → matrix → 04_vmr_scanpy.py`：对 0.01/0.02/0.05 三个 VMR threshold 分支分别运行。
10. `05_summarize_run.py`: validates pairwise meth-diff, the hypo-DMR figures, and the three VMR branches, then generates the summary. / `05_summarize_run.py`：验证 pairwise meth-diff、hypo-DMR 图和三个 VMR 分支后生成 summary。

## Output structure / 输出结构

```text
<run>/
├── 00_scanpy_selected/
├── 01_cov/cov/*.cov.gz
├── 02_prepared/
├── 03_filtered/
├── 07_methdiff/samples/<sample>/comparisons/<cell_type_A>_vs_<cell_type_B>/DMRs.bed
├── 07_methdiff/pairwise_summary.tsv
├── 07_methdiff/pairwise_summary.json
├── 08_rawp_fallback/fallback_summary.json
├── 09_top200_hypo_DMRs/samples/<sample>/
├── 10_single_cell_DMR_matrix/samples/<sample>/single_cell_DMR_mean_unique_CpG_ratio.tsv.gz
├── 11_hypo_DMR_heatmaps/samples/<sample>/mean_ratio.png
├── 11_hypo_DMR_heatmaps/samples/<sample>/zscore.png
├── 11_hypo_DMR_heatmaps/samples/<sample>/zscore_colorbar_compressed.png
├── pooled/07_methdiff/samples/<pooled-label>/comparisons/
├── pooled/07_methdiff/pairwise_summary.tsv
├── pooled/07_methdiff/pairwise_summary.json
├── 04_scan/var_0.01|var_0.02|var_0.05/VMRs.bed
├── 05_matrix/var_0.01|var_0.02|var_0.05/
├── 06_scanpy/var_0.01|var_0.02|var_0.05/
├── stage_status.tsv
├── <stage>.resources.json
├── run_summary.json
└── run_summary.tsv
```

`pooled/07_methdiff/samples/<pooled-label>/` uses the pooled label as the sample directory name (example/示例: `CYL_ZCP`). / `pooled/07_methdiff/samples/<pooled-label>/` 以 pooled 标签作为样本目录名（示例/example: `CYL_ZCP`）。

Main-pipeline completion criteria: `run_summary.json` has `status=complete`; the per-sample pairwise/fallback/Top200/matrix/heatmap summaries are complete; the h5ad, UMAP, and matrix of all three VMR branches exist and the cell-ID check passes. The pooled branch is judged independently by `pooled/methdiff.COMPLETE` and `status=complete` in `pooled/07_methdiff/pairwise_summary.json`. Cancelled or failed directories are retained for diagnosis and are not rerun in place by overwriting. / 主流程完成判据：`run_summary.json` 的 `status=complete`，分样本 pairwise/fallback/Top200/matrix/heatmap summary 完成，三个 VMR 分支的 h5ad、UMAP、matrix 均存在且 cell ID 校验通过。pooled 分支独立以 `pooled/methdiff.COMPLETE` 和 `pooled/07_methdiff/pairwise_summary.json` 的 `status=complete` 为准。取消或失败的目录保留用于诊断，不覆盖后重跑。
