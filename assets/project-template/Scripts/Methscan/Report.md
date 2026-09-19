# MethSCAn workflow report / MethSCAn 流程报告（整理版）

## Script numbering / 脚本编号

The previous version had `02_prepare_methscan.py` and `02_convert_allc_to_cov.py` sharing the same number; they were renumbered by actual dependency order: / 旧版 `02_prepare_methscan.py` 与 `02_convert_allc_to_cov.py` 重号，已按实际依赖重排为：

```text
01_select_scanpy_cells.py
02_convert_allc_to_cov.py
03_prepare_methscan.py
04_vmr_scanpy.py
05_summarize_run.py
06_run_with_resources.py
07_methdiff_celltype.py
08_hypo_dmr_heatmaps.py
```

The main entry `run_methscan.sbatch` was updated in step: it no longer calls the non-existent `01_prepare_allc_inputs.py` and no longer uses the old `link_path`/`00_manifest` intake stages. / 主入口 `run_methscan.sbatch` 已同步更新引用，不再调用不存在的 `01_prepare_allc_inputs.py`，也不再使用旧的 `link_path`/`00_manifest` intake 阶段。

Added: `run_methscan_common.sbatch`, `run_methscan_qc_stage.sbatch`, `run_methscan_branch.sbatch`, `run_methscan_summary.sbatch`, and `submit_methscan_pipeline.sh`. A later production run then needs only one `submit_methscan_pipeline.sh` call: it automatically submits the common stages, three independent QC stages, three threshold branches, and the summary, and sets the `afterok` dependencies. prepare/filter/smooth each request 4 CPU/16G. / 新增 `run_methscan_common.sbatch`、`run_methscan_qc_stage.sbatch`、`run_methscan_branch.sbatch`、`run_methscan_summary.sbatch` 和 `submit_methscan_pipeline.sh`。后续正式运行只需执行一次 `submit_methscan_pipeline.sh`：它自动提交公共阶段、三个独立 QC 阶段、三个 threshold 分支和 summary，并设置 `afterok` 前后依赖。prepare/filter/smooth 各申请 4 CPU/16G。

## Current workflow / 当前流程

```text
ALLC → Scanpy whitelist → ALLC→cov → prepare → filter → smooth
→ per-sample pairwise cell-type meth-diff          / 样本内 pairwise cell-type meth-diff
→ optional raw-p fallback                          / 可选 raw-p fallback
→ Top200 hypo-DMR                                  / Top200 hypo-DMR
→ single-cell × DMR matrix                          / 单细胞×DMR矩阵
→ raw mean ratio, z-score, compressed-colorbar z-score heatmaps / raw mean ratio / z-score / 压缩 colorbar z-score 热图
→ (scan → matrix → Scanpy) × {0.01, 0.02, 0.05} → summary

independent optional branch / 独立可选分支:
smooth → merge same-named cell types across samples / 合并各样本同名 cell type
→ <N> pooled pairwise meth-diff comparisons (example/示例: 91) / <N> 个 pooled pairwise meth-diff（示例/example: 91）
```

Pairwise meth-diff uses a hard-link view of `03_filtered` after smoothing, restricted to the chromosomes declared by `references.chrom_sizes` (no data is copied; the example/示例 project declares chr1–22/X/Y and the summary records the scope actually used). The main-pipeline group files restrict cells per sample and write results under `07_methdiff/samples/<sample>/comparisons/`; the independent pooled branch merges same-named cell types across the samples (example/示例: `CYL`, `ZCP`) and writes under `pooled/07_methdiff/samples/<pooled-label>/` (example/示例: `CYL_ZCP`). The per-sample heatmap downstream strictly follows the reference workflow's convention of raw p `<0.01`, absolute mean-ratio difference `>=0.25`, and Top200 hypo-DMRs per cell type; the matrix comes from the equal-weight mean of the unique CpG ratio inside each COV. The three heatmap sets are the raw mean ratio, the DMR-wise z-score, and the value-preserving z-score whose colorbar alone is compressed to `[-1,1]`; they display row/column cell-type labels and group boundaries. / pairwise meth-diff 使用 smooth 后 `03_filtered` 的 chr1–22/X/Y hard-link view（不复制数据）。主流程 group 文件按 sample 限定细胞，在 `07_methdiff/samples/<sample>/comparisons/` 下写入结果；独立 pooled 分支把各样本同名 cell type 合并后写入 `pooled/07_methdiff/samples/<pooled-label>/`（示例/example: `CYL_ZCP`）。分样本热图下游严格采用参考流程的 raw p `<0.01`、mean ratio 绝对差 `>=0.25`、每细胞类型 Top200 hypo-DMR 口径；矩阵来自 COV 内 unique CpG ratio 的等权平均。三套热图分别为原始 mean ratio、DMR-wise z-score，以及仅把 colorbar 压缩到 `[-1,1]` 的未改值 z-score，并显示行列 cell-type 标签与分组边界。

The ALLC→cov stage converts cells with multiple processes; `run_methscan.sbatch` sets the worker count to the Slurm-allocated CPU count, and `SCMO_COV_WORKERS` can override it manually. The current production job had already completed conversion with 16 workers before this change; later jobs automatically follow the requested CPU count. / ALLC→cov 阶段按细胞使用多进程转换；`run_methscan.sbatch` 将 worker 数设置为 Slurm 分配的 CPU 数，也可通过 `SCMO_COV_WORKERS` 手动覆盖。当前正式任务在该修改前已使用 16 workers 完成转换；后续新任务会自动跟随申请的 CPU 数。

## Key filtering rules / 关键筛选规则

- canonical cell ID: `<sample_id>_<17bp_barcode>`; / canonical cell ID：`<sample_id>_<17bp_barcode>`；
- `cell_type` empty, blank, or `NA`: excluded; / `cell_type` 为空、空白或 `NA`：排除；
- ALLC not present in the RNA annotation table: excluded; / 不在 RNA 注释表中的 ALLC：排除；
- the ALLC source directory is read-only; conversion results are written into the current run's `01_cov`. / ALLC 源目录只读，转换结果写入当前 run 的 `01_cov`。

## Run and verification / 运行与验证

Submission follows the scheduler profile's preferred node; the resource numbers quoted in this section are an illustrative profile, and the actual values come from the scheduler profile. / 提交按 scheduler profile 的首选节点；本节引用的资源数字为示例 profile，实际由 scheduler profile 决定：

```bash
sbatch <scheduler-profile-resource-options> \
  "$SCMO_PROJECT_ROOT/Scripts/Methscan/run_methscan.sbatch" \
  "$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>"
```

Split submission entry point / 拆分提交入口：

```bash
"$SCMO_PROJECT_ROOT/Scripts/Methscan/submit_methscan_pipeline.sh" \
  "$SCMO_PROJECT_ROOT/Results/Methscan/<run-name>"
```

The common job requests 55 CPU/250G by default; each threshold branch requests 18 CPU/80G by default and is scheduled by Slurm onto an available CPU node. / 公共作业默认申请 55 CPU/250G；每个 threshold 分支默认申请 18 CPU/80G，由 Slurm 自动调度到可用 CPU 节点。

prepare, filter, and smooth have been split into independent serial jobs, each requesting 4 CPU/16G. Pairwise meth-diff runs rolling at `2 comparisons × 4 threads` by default and requests 8 CPU/16G; `07_methdiff_celltype.py` enforces `jobs × threads ≤ SLURM_CPUS_PER_TASK`, and because each comparison uses its own process, its own output, and MethSCAn's fixed random seed, concurrency does not change the statistical convention of a single comparison. The hypo-DMR matrix/plotting stage uses 8 workers, 8 CPU/48G by default. Both may be scheduled into whichever of the profile's partitions has free capacity (example/示例: `cpu`, `fat`), and can be overridden through `SCMO_METHDIFF_JOBS`, `SCMO_METHDIFF_THREADS`, `SCMO_HYPO_MATRIX_WORKERS`, or `sbatch` resource options. selection+ALLC→COV uses the illustrative 55 CPU/250G profile with an unpinned node. / prepare、filter、smooth 已拆为独立串行作业，各申请 4 CPU/16G。pairwise meth-diff 默认以 `2 comparisons × 4 threads` 滚动执行，申请 8 CPU/16G；`07_methdiff_celltype.py` 强制检查 `jobs × threads ≤ SLURM_CPUS_PER_TASK`，各比较使用独立进程、独立输出且 MethSCAn 固定随机种子，因而并发不改变单比较统计口径。hypo-DMR 矩阵/绘图默认使用 8 workers、8 CPU/48G。两者允许在 profile 的分区（示例/example: `cpu,fat`）间择空调度，并可通过 `SCMO_METHDIFF_JOBS`、`SCMO_METHDIFF_THREADS`、`SCMO_HYPO_MATRIX_WORKERS` 或 `sbatch` 资源参数覆盖。selection+ALLC→COV 使用示例 profile 的 55 CPU/250G，节点不固定。

Check `stage_status.tsv`, the Slurm logs, and each `<stage>.resources.json`. The following must hold at the end: / 检查 `stage_status.tsv`、Slurm 日志和每个 `<stage>.resources.json`。最终必须满足：

- `run_summary.json` exists and `status=complete`; / `run_summary.json` 存在且 `status=complete`；
- `07_methdiff/pairwise_summary.json` contains no hard failures; comparisons with the known FDR division-by-zero are rerun in a targeted way by `08_rawp_fallback`; / `07_methdiff/pairwise_summary.json` 不含硬失败；已知 FDR 除零比较由 `08_rawp_fallback` 定向补跑；
- `hypo_dmr_heatmap_summary.json` exists and `status=complete`; / `hypo_dmr_heatmap_summary.json` 存在且 `status=complete`；
- all three threshold branches have VMRs, the four matrices, Scanpy h5ad, and figures; / 三个 threshold 分支均有 VMR、四个 matrix、Scanpy h5ad 和图；
- matrix row IDs match the filtered header; / matrix 行 ID 与 filtered header 一致；
- the Scanpy embedding is a subset of the filtered cells; / Scanpy embedding 是 filtered cell 的子集；
- all resource records have `return_code=0`. / 所有资源记录 `return_code=0`。

Partially run or cancelled directories must not be resumed in place; after repair, use a new output directory. / 部分运行或取消的目录不得直接续跑；修复后使用新的输出目录。

## Production run record / 正式运行记录

Job `<job-id>` completed the whole `<run-name>` pipeline on `<date>`. Of `<N>` ALLC files (example/示例: 8,949), `<N>` passed the RNA cell-ID and non-empty/non-`NA` cell type filter (example/示例: 8,626); MethSCAn filter retained `<N>` cells (example/示例: 6,264). The three branches produced `<N>`, `<N>`, and `<N>` VMRs respectively (example/示例: 39,553, 80,818, 166,618); all three branches completed their matrix and Scanpy outputs and `run_summary.json` was `status=complete`. prepare finished at `<date> <time>` and the final summary at `<date> <time>`. / 作业 `<job-id>` 于 `<date>` 完成 `<run-name>` 全流程。`<N>` 个 ALLC 中 `<N>` 个通过 RNA cell-ID 及非空/非 `NA` cell type 筛选（示例/example: 8,949 中 8,626 个）；MethSCAn filter 后保留 `<N>` 个细胞（示例/example: 6,264）。三个分支分别产生 `<N>`、`<N>`、`<N>` 个 VMR（示例/example: 39,553、80,818、166,618），三个分支均完成 matrix 和 Scanpy 输出，`run_summary.json` 为 `status=complete`。prepare 于 `<date> <time>` 完成，最终 summary 于 `<date> <time>` 完成。

The per-sample hypo-DMR downstream has completed: sample A produced `<N>` merged Top200 hypo-DMRs × `<N>` cells (example/示例: 2,409 × 2,919) and sample B produced `<N>` × `<N>` (example/示例: 2,502 × 3,345); both generated mean-ratio, z-score, and compressed-colorbar z-score heatmaps with cell-type row/column labels. / 分样本 hypo-DMR 下游已完成：样本 A 为 `<N>` 个合并 Top200 hypo-DMR × `<N>` 个细胞（示例/example: 2,409 × 2,919），样本 B 为 `<N>` × `<N>`（示例/example: 2,502 × 3,345）；两者均生成带 cell-type 行列标签的 mean-ratio、z-score 和压缩 colorbar z-score 热图。

Pooled job `<job-id>` completed on `<date>` after `<duration>`. It merged same-named RNA cell types across the samples and completed `<N>`/`<N>` pairwise comparisons among `<N>` cell types over all `<N>` filtered cells (example/示例: 91/91 comparisons, 14 cell types, 6,264 cells); `<N>` failures, `<N>` fallbacks, `<N>` empty results (example/示例: 0, 0, 0), and `<N>` DMR rows in total (example/示例: 6,437,767). The pooled stage is responsible for MethSCAn DMRs only and generates no Top200, single-cell matrix, or heatmaps; its full unique-DMR output is later consumed by `Scripts/Methylvi/vmr_dmr/`. / pooled 作业 `<job-id>` 于 `<date>` 完成，耗时 `<duration>`。它将各样本的同名 RNA cell type 合并，在全部 `<N>` 个过滤后细胞中完成 `<N>` 种细胞类型的 `<N>`/`<N>` 个两两比较（示例/example: 6,264 个细胞中 14 种细胞类型的 91/91 个比较）；`<N>` 失败、`<N>` fallback、`<N>` 空结果（示例/example: 0、0、0），共 `<N>` 行 DMR（示例/example: 6,437,767）。pooled 阶段只负责 MethSCAn DMR，不生成 Top200、单细胞矩阵或热图；其全量 unique-DMR 后续由 `Scripts/Methylvi/vmr_dmr/` 消费。
