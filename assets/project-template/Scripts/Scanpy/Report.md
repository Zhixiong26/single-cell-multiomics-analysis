# Scanpy notebook iteration report / Scanpy Notebook 迭代报告

This file is a continuously updated run ledger. The latest state comes first; historical records keep their conclusions and evidence and do not copy the notebook's full implementation. Stable rules live in `README.md`.

本文件是可持续更新的运行账本。最新状态放在最前；历史记录保留结论和证据，不复制 Notebook 中的完整实现。稳定规则见 `README.md`。

## 最新状态 / Latest status

- 状态 / Status：`<运行后填写 / record at run time>`（示例 / example: 已实现优先改进；baseline 与隔离 candidate smoke test 均已无错误完成 / the priority improvements are implemented; baseline and the isolated candidate smoke test both completed without errors）。
- Canonical Notebook：`Notebooks/<notebook>.ipynb`（打包入口 / packaged entry: `Notebooks/scanpy_workflow.ipynb`）。
- 本轮 Notebook 行为 / Notebook behaviour this round：`<填写 / fill in>`（示例 / example: 已修改并完整执行；17 个代码单元顺序执行且无 error output / modified and executed in full; 17 code cells ran in order with no error output）。完整执行副本保存到 `Results/Scanpy/<run-id>/`，canonical Notebook 随后清为零输出交付状态。
  The complete executed copy is saved into the results directory and the canonical notebook is then cleared back to its zero-output delivery state.
- 正式输出目录 / formal output directory：`Results/Scanpy/<run-id>/`（`<run-id>` 为本轮运行标识，例如 `<sample-pair>_notebook` / `<run-id>` identifies this run, e.g. `<sample-pair>_notebook`）。
- 已知保存开关 / known save switches：`ANALYSIS_CONFIRMED`、`OVERWRITE_DATA_OUTPUTS`（模板中由运行参数 sidecar `$SCMO_SCANPY_PARAMETERS` 注入，源自 `config/analysis.yaml` 的 `analysis.scanpy.notebook`）；每次运行前必须重新读取确认。示例值 / example values: `ANALYSIS_CONFIRMED = True`、`OVERWRITE_DATA_OUTPUTS = False`。
  In the template these are injected by the run parameter sidecar `$SCMO_SCANPY_PARAMETERS`, sourced from `analysis.scanpy.notebook` in `config/analysis.yaml`; the values must be re-read and confirmed before every run.
- 已授权范围 / authorized scope：可修改并重跑 PCA、Harmony、邻接图、UMAP、Leiden 参数；候选运行不得覆盖正式结果。
  The PCA, Harmony, neighbour-graph, UMAP and Leiden parameters may be modified and re-run; a candidate run must not overwrite formal results.
- 当前未决事项 / open items：暂无实现阻塞；真实参数优化时应使用新的 `ITERATION_ID`，并比较候选 UMAP、marker 和审计表。
  No implementation blocker at present; a real parameter optimization must use a new `ITERATION_ID` and compare the candidate UMAP, markers and audit tables.

## 已知参考基线 / Known reference baseline

The following are replay records of the existing notebook, used for difference detection; they are not mandatory success conditions.

以下是既有 Notebook 重放记录，用于差异检测，不是强制成功条件：

| 项目 / Item | 已记录值 / Recorded value | 使用方式 / How it is used |
|---|---:|---|
| CYL 输入 cells / input cells | `<N>`（示例 / example: 4,264） | 输入一致时比较 / compare when the input is unchanged |
| ZCP 输入 cells / input cells | `<N>`（示例 / example: 4,685） | 输入一致时比较 / compare when the input is unchanged |
| QC 与 doublet 后 cells / cells after QC and doublet removal | `<N>`（示例 / example: 8,696） | 差异出现时定位 QC/Scrublet 阶段 / locate the QC/Scrublet stage when a difference appears |
| Leiden clusters | `<N>`（示例 / example: 18（0–17）） | 仅作结构变化警报；不得据此强套标签 / structural change alarm only; labels must never be forced from this number |
| 随机种子 / random seed | 0 | 与实际参数 JSON 交叉检查 / cross-check against the actual parameter JSON |

Recorded formal artifacts (`CYL`/`ZCP` in the example name are sample labels, 示例 / example):

已记录的正式产物：

- `Results/Scanpy/<run-id>/<output_stem>_confirmed.h5ad`（示例 / example: `rna_e_cyl_zcp_confirmed.h5ad`）
- `Results/Scanpy/<run-id>/cell_id_cell_type.tsv`
- `Results/Scanpy/<run-id>/confirmed_parameters.json`
- `Results/Scanpy/<run-id>/figure_manifest.json`
- `Results/Scanpy/<run-id>/figures/`

These paths express the expected interface; they do not prove that the files still exist or are still valid at the next run.

这些路径表示预期接口，不证明文件在下一次运行时仍存在或仍有效。

## 当前验证矩阵 / Current validation matrix

| 层级 / Layer | 状态 / Status | 当前证据 / Current evidence | 下一步 / Next step |
|---|---|---|---|
| 执行完整性 / Execution integrity | `<未验证 / 通过 / 失败 / 阻塞>`（示例 / example: 通过） | `<本轮证据 / this round's evidence>`（示例 / example: baseline 与 candidate 的 17 个代码单元均按 1–17 执行，无 error output） | 真实候选继续同样检查 / apply the same check to real candidates |
| 结构完整性 / Structural integrity | `<未验证 / 通过 / 失败 / 阻塞>` | `<本轮证据>`（示例 / example: baseline 与 candidate 均为 8,696 cells、18 clusters；candidate H5AD/TSV 结构一致） | 参数改变后重新检查 / re-check after any parameter change |
| 分析合理性 / Analytical plausibility | `<未验证 / 通过 / 失败 / 阻塞>` | `<本轮证据>`（示例 / example: baseline 重放生成 20 张图；核心细胞数和 cluster 结构保持一致） | 后续候选需单独比较 UMAP/marker / later candidates compare UMAP and markers individually |
| 注释完整性 / Annotation completeness | `<未验证 / 通过 / 失败 / 阻塞>` | `<本轮证据>`（示例 / example: baseline 为 `reviewed`；candidate 生成 18 行审计，17 个 proposed、1 个 ambiguous/Unassigned，状态为 `candidate_requires_review`） | 候选不得直接升级为正式标签 / a candidate must not be promoted to formal labels directly |
| 输出一致性 / Output consistency | `<未验证 / 通过 / 失败 / 阻塞>` | `<本轮证据>`（示例 / example: baseline manifest 20 张图；candidate manifest 18 张图且零缺失；正式三个数据文件哈希在 smoke test 前后不变） | 每个真实候选使用新目录 / use a new directory for every real candidate |

状态只使用：`未验证`、`通过`、`失败`、`阻塞`。没有证据时不得写为“通过”。

Only these statuses may be used: `未验证` (unverified), `通过` (pass), `失败` (fail), `阻塞` (blocked). Without evidence, nothing may be written as "pass".

## 迭代记录 / Iteration log

New entries go at the top of this section, using the template in "新一轮记录模板". The three entries below are example records carried over from the porting source: their dates, job evidence and measured numbers are replaced by placeholders that show what belongs there; a new project replaces them with its own records.

新条目写在“迭代记录”顶部，使用“新一轮记录模板”。以下三条为移植来源的示例记录：日期、作业证据和测量数字已替换为占位符，以显示该位置应填写什么；新项目应替换为本轮真实记录。

### YYYY-MM-DD（示例 / example: 2026-08-27）：实现候选参数与注释自迭代后端 / implement the candidate-parameter and annotation self-iteration backend

- 触发 / Trigger：终检发现错误输出污染、固定输出目录、UMAP 参数未集中和 cluster 变化后缺少候选注释流程。
  The final check found polluted error outputs, a fixed output directory, UMAP parameters that were not centralized, and no candidate annotation flow after a cluster change.
- 修改 / Change：新增 baseline/candidate 运行模式、唯一迭代目录、集中式 PCA/Harmony/neighbors/UMAP/Leiden 参数、marker-score 候选注释与审计输出。
  Added baseline/candidate run modes, unique iteration directories, centralized PCA/Harmony/neighbors/UMAP/Leiden parameters, and marker-score candidate annotation with audit output.
- Notebook 清洁 / Notebook cleanup：清除 `<N>` 个 `notebook controller is DISPOSED` 输出（示例 / example: 16），并用 `<运行环境 / run environment>` kernel 完整重放（模板中为 `environments.tsv` 声明的 `scanpy_allcools` 解释器 / in the template, the `scanpy_allcools` interpreter declared in `environments.tsv`）。
- baseline 验证 / baseline verification：`<N>` 个代码单元顺序完成、无错误（示例 / example: 17）；三个正式数据文件通过逐细胞保护而未覆盖；`<N>` 张图重新生成（示例 / example: 20）。
- 执行证据 / Execution evidence：完整副本保存为 `Results/Scanpy/<run-id>/scanpy_workflow_executed.ipynb`；canonical Notebook 验证后再次清除输出，避免提交内嵌图和 kernel 状态。
- 结果身份 / Result identity：manifest 写入 `run_kind=baseline`、`iteration_id=null`、`annotation_status=reviewed`。
- candidate smoke test：使用相同分析参数但强制进入 candidate 分支，独立写入 `<独立候选目录 / independent candidate directory>`（示例 / example: `/tmp/<candidate_smoke>/iterations/same_params_smoke/`）。
- candidate 结果 / candidate result：`<N>` 个 cluster 得到 marker-score proposal（示例 / example: 17）；cluster `<id>` 因 `<top1>/<top2>` top1-top2 margin 仅 `<margin>`（示例 / example: cluster 13, Mast/T, 0.058），小于 0.20 阈值而保守标为 `Unassigned`。
- 隔离验证 / Isolation check：candidate 使用 `candidate_cell_type`，没有正式 `cell_type` 列；H5AD/TSV/JSON 和 `<N>` 张图均使用候选文件名（示例 / example: 18），正式数据哈希未变化。
- 结论 / Conclusion：baseline 保护、候选隔离、候选审计与低置信拒绝机制均工作正常。
  baseline protection, candidate isolation, candidate auditing and the low-confidence rejection mechanism all work as intended.

### YYYY-MM-DD（示例 / example: 2026-08-27）：扩展注释 UMAP 自迭代范围 / extend the self-iteration scope to the annotation UMAP

- 触发 / Trigger：用户要求迭代不能只处理报错和重复确认注释，还要评价最终注释 UMAP，并在存在明显优化空间时自主调参、重跑和再次判断。
  The user required that iteration not merely handle errors and re-confirm annotations, but also evaluate the final annotation UMAP and autonomously tune parameters, re-run and judge again when an obvious opportunity exists.
- 授权 / Authorization：允许修改 PCA、Harmony、邻接图、UMAP 和 Leiden 参数；不包含静默改变输入、QC、生物学结论或覆盖正式结果。
  Modification of the PCA, Harmony, neighbour-graph, UMAP and Leiden parameters is allowed; it does not include silently changing inputs, QC, biological conclusions, or overwriting formal results.
- 修改 / Change：README 新增 UMAP 合理性检查、可调参数家族、受影响重跑范围、候选比较、接受标准和收敛条件。
  Added to the README: the UMAP plausibility checks, the adjustable parameter families, the affected re-run scope, candidate comparison, acceptance criteria and convergence conditions.
- 验证 / Verification：本轮仅更新文档，没有执行 Notebook 或产生候选分析。
  This round only updated documentation; the notebook was not executed and no candidate analysis was produced.
- 结论 / Conclusion：后续 skill 应持续迭代到分析证据收敛，而不是以无报错或已有 cell-type 标签作为终点。
  The skill must keep iterating until the analytical evidence converges, rather than treating "no error" or an existing cell-type label as the end point.

### YYYY-MM-DD（示例 / example: 2026-08-27）：文档结构修正 / documentation structure correction

- 触发 / Trigger：需要让 README/Report 可被 skill 稳定使用，并支持自我迭代和纠错。
  The README and Report had to be usable reliably by a skill and to support self-iteration and error correction.
- 观察 / Observation：旧文档能说明入口和输出，但没有失败分类、修正边界、停止条件或统一记录格式。
  The old documentation described the entry point and outputs, but had no failure classification, correction boundaries, stop conditions or unified record format.
- 修改 / Change：README 改为稳定操作契约；Report 改为状态账本，增加参考基线、验证矩阵、修正记录和迭代模板。
  The README became a stable operational contract; the Report became a status ledger with a reference baseline, validation matrix, correction records and an iteration template.
- 验证 / Verification：核对 Notebook 中的输出目录、随机种子、保存开关和正式输出文件名；未修改 Notebook。
  Checked the output directory, random seed, save switches and formal output filenames in the notebook; the notebook was not modified.
- 结论 / Conclusion：文档层面的缺口已修正；分析结果本轮未重新验证。
  The documentation gap is corrected; the analysis results were not re-verified this round.

## 修正决策账本 / Correction decision ledger

| 日期 / Date | 现象 / Symptom | 根因证据 / Root-cause evidence | 最小修正 / Minimal correction | 再验证结果 / Re-verification result | 是否需人工决定 / Needs a human decision |
|---|---|---|---|---|---|
| YYYY-MM-DD（示例 / example: 2026-08-27） | 迭代未覆盖最终 UMAP 优化 / iteration did not cover final UMAP optimization | 用户明确要求审查并自主调整降维/聚类参数 / the user explicitly required reviewing and autonomously adjusting the dimensionality-reduction/clustering parameters | 扩展 README 迭代规则与权限 / extend the README iteration rules and authorization | Markdown 与规则一致性检查通过；分析未运行 / markdown and rule-consistency checks passed; the analysis was not run | 否 / no |
| YYYY-MM-DD（示例 / example: 2026-08-27） | 文档不支持闭环迭代 / the documentation did not support closed-loop iteration | 缺少验证层级、停止条件和记录模板 / the validation layers, stop conditions and record template were missing | 仅重构 README/Report / restructure README/Report only | 文档检查通过；分析未运行 / documentation checks passed; the analysis was not run | 否 / no |

## 新一轮记录模板 / New-round record template

Copy the section below to the top of "迭代记录" (Iteration log) and fill in the actual evidence; delete the inapplicable items and do not leave placeholder sentences behind.

复制以下小节到“迭代记录”顶部，并填写实际证据；删除不适用项，不保留占位语句。

```markdown
### YYYY-MM-DD：<本轮目标 / round objective>

- 请求与权限 / Request and permission：<允许读取、执行、修改或覆盖的范围 / what may be read, executed, modified or overwritten>
- Notebook SHA-256 / Git 状态：<值 / value>
- 环境 / Environment：<kernel、Python、Scanpy、关键依赖 / kernel, Python, Scanpy, key dependencies>
- 输入 / Inputs：<路径、身份、维度、校验信息 / paths, identity, dimensions, checksums>
- 现有输出 / Existing outputs：<文件及是否允许覆盖 / files and whether overwriting is allowed>
- 执行结果 / Execution result：<完成位置、错误、耗时和资源信息 / completion point, errors, duration and resource information>
- 五层验证 / Five-layer verification：<执行、结构、分析、注释、输出 / execution, structure, analysis, annotation, output>
- 注释判定 / Annotation calls：<逐 cluster 的 cell type 与依据 marker，以及标记为不确定的 cluster / the cell type adopted per cluster with the markers behind it, and which clusters are flagged uncertain>
- 注释确认 / Annotation confirmation：<用户已确认 / 已请求确认 / 待修正；修正内容与重跑范围 / confirmed by the user, asked and pending, or corrected — with what changed and what was re-run>
- 最终注释 UMAP 诊断 / Final annotation UMAP diagnosis：<碎裂、混合、样本效应、QC 驱动、稀有群、marker/邻域一致性 / fragmentation, mixing, sample effects, QC-driven structure, rare populations, marker/neighbourhood agreement>
- 参数候选与假设 / Parameter candidates and hypotheses：<本轮只改变的参数家族、候选值、预期改善 / the single parameter family changed this round, candidate values, expected improvement>
- 候选输出目录 / Candidate output directory：<独立且不覆盖正式结果的路径 / an independent path that does not overwrite formal results>
- 与上一基线比较 / Comparison with the previous baseline：<改善、退化、稳定性和是否接受 / improvement, degradation, stability and whether it was accepted>
- 与参考基线差异 / Differences from the reference baseline：<差异及其解释，不以匹配为目标 / the difference and its explanation; matching is not the goal>
- 根因证据 / Root-cause evidence：<失败时填写 / fill in on failure>
- 最小修正 / Minimal correction：<修改内容；未修改则写明 / what was changed; state explicitly if nothing was>
- 再验证 / Re-verification：<从何处重跑、哪些下游检查通过 / where the re-run started and which downstream checks passed>
- 结论 / Conclusion：<通过、失败或阻塞 / pass, fail or blocked>
- 未决风险与下一步 / Open risks and next step：<需要用户决定的事项 / items that need the user's decision>
```

## 更新规则 / Update rules

- New state is written at the top; historical conclusions are never silently rewritten.
  新状态写在顶部，不静默改写历史结论。
- Failure evidence, anomalous tracebacks and difference summaries must keep locatable paths.
  失败证据、异常 traceback 和差异摘要应保留可定位路径。
- Parameter values are governed by the notebook and `confirmed_parameters.json`; the report records only this round's key differences.
  参数值以 Notebook 与 `confirmed_parameters.json` 为准；报告只记录本轮关键差异。
- After a change, check markdown, paths, factual consistency and the Git diff.
  修改完成后检查 Markdown、路径、事实一致性和 Git diff。
- If new evidence overturns an old conclusion, add a correction record stating what was corrected, the evidence and the affected scope.
  若新证据推翻旧结论，新增一条更正记录，说明被更正内容、证据和影响范围。
