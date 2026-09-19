# Environment 分析报告 / Analysis report

状态 / Status：`<待填写 / fill in>`（示例/example: 环境发现与初始映射已完成 / discovery and initial mapping complete）

## 已验证环境 / Verified environments

每个已验证环境记录：环境名、环境根目录（`$SCMO_ENV_ROOTS`）、关键工具及版本，以及用于验证的 profile 验证命令。 / For each verified environment record the environment name, its root under `$SCMO_ENV_ROOTS`, the key tools with versions, and the profile verification command that was used.

- BAM/ALLCools/Scanpy：`<环境名 / env name>`，Samtools `<版本 / version>`、ALLCools `<版本 / version>`、Scanpy `<版本 / version>`。
  - 示例 / example：`allcools`，Samtools 1.22、ALLCools 1.1.1、Scanpy 1.9.3。
- Methscan：`<环境名 / env name>`，MethSCAn `<版本 / version>`、Python `<版本 / version>`。
  - 示例 / example：`MethSCAn`，MethSCAn 1.1.0、Python 3.8.20。
- MethylVI：`<环境名 / env name>`，Scanpy `<版本 / version>`、scvi-tools `<版本 / version>`、Torch `<版本 / version>`。
  - 示例 / example：`methylvi`，Scanpy 1.12.3、scvi-tools 1.5.0.post1、Torch 2.13.0。
- FASTQ/Bismark 候选 / candidate：`<环境名 / env name>`，Bismark `<版本 / version>`、Bowtie2 `<版本 / version>`，以及不可用于该用途的工具。
  - 示例 / example：`moabs`，Bismark 0.24.2、Bowtie2 2.5.4；其中 Samtools 0.1.19 过旧，不用于 BAM 验证。 / Samtools 0.1.19 in that environment is too old and is not used for BAM verification.

环境根目录 / Environment roots：`$SCMO_ENV_ROOTS`（示例/example: `$HOME/miniconda3/envs`、`$HOME/miniforge3/envs`，并可通过 `CONDA_ENVS_PATH` 追加）。

所选混合映射记录在 `config/environments.tsv`。本项目之外的共享环境是只读依赖，不得修改。 / The selected hybrid mapping is recorded in `config/environments.tsv`. Shared environments outside this project are read-only dependencies and must not be modified.

## 发现结果 / Discovery results

- 完整只读清单位于 `Results/Environment/environment_inventory.tsv`。 / The complete read-only inventory is at `Results/Environment/environment_inventory.tsv`.
- 初次发现脚本因旧版 Bash 对空数组的 `set -u` 行为失败，已修正并重新完整运行。 / The first discovery attempt failed because of old-Bash empty-array behavior under `set -u`; it was fixed and rerun successfully.
- 未创建、安装、升级或删除任何环境。 / No environment was created, installed into, upgraded, or removed.
- `<其他发现与失败 / other findings and failures>`

## 服务器资源快照 / Server resource snapshot

检查时间 / Checked：`<日期与时间 / date and time>`（示例/example: 2026-08-24 14:23 CST）

- 登录节点 `<登录节点名 / login node>`：`<N>` 个物理 CPU cores（示例/example: 32；2 × Xeon Silver 4314），`<N>` GiB RAM，其中约 `<N>` GiB available（示例/example: 125 GiB，约 104 GiB available）；swap `<N>` GiB，已使用约 `<N>` GiB（示例/example: 31 GiB，约 15 GiB）。登录节点仅用于轻量检查，不作为完整 QC 运行资源。 / Login node `<login node>`: `<N>` physical CPU cores (example: 32; 2 × Xeon Silver 4314), `<N>` GiB RAM with about `<N>` GiB available (example: 125 GiB, about 104 GiB available); `<N>` GiB swap with about `<N>` GiB used (example: 31 GiB, about 15 GiB). The login node is for lightweight checks, not full QC runs.
- 项目文件系统 `$SCMO_PROJECT_ROOT`：`<N>` TiB 中仅约 `<N>` GiB available，使用率 `<N>`%（示例/example: 204 TiB 中约 1 GiB available，100%）；inode 使用率 `<N>`%（示例/example: 96）。在清理或确定新结果位置前，不应向这里写入大型中间文件。 / Project filesystem `$SCMO_PROJECT_ROOT`: only about `<N>` GiB available out of `<N>` TiB, `<N>`% used (example: about 1 GiB out of 204 TiB, 100%); inode usage is `<N>`% (example: 96). Do not write large intermediates here until space is cleared or a new output location is selected.
- 源数据盘 `<源数据文件系统 / source-data filesystem>`：`<N>` TiB，总体使用率 `<N>`%，约 `<N>` TiB available（示例/example: 100 TiB，96%，约 4.7 TiB）；inode 余量充足。 / Source-data filesystem `<source-data filesystem>`: `<N>` TiB total, `<N>`% used, about `<N>` TiB available (example: 100 TiB, 96%, about 4.7 TiB); inode capacity is ample.
- Slurm `cpu`：`<N>` × `<N>`-core nodes、每节点约 `<N>` GiB RAM（示例/example: 3 × 56-core nodes，约 257.4 GiB）；快照时 `<N>`/`<N>` CPUs allocated（示例/example: 23/168）。`<节点名>` 与 `<节点名>` 为 idle（示例/example: `cu02`、`cu03`）。 / Slurm `cpu`: `<N>` × `<N>`-core nodes with about `<N>` GiB RAM each (example: three 56-core nodes, about 257.4 GiB each); `<N>`/`<N>` CPUs were allocated at the snapshot (example: 23/168). Nodes `<node-a>` and `<node-b>` were idle (example: `cu02`, `cu03`).
- Slurm `fat`：`<N>` × `<N>`-core node、约 `<N>` RAM（示例/example: 1 × 192-core node，约 1.03 TiB）；快照时 `<N>`/`<N>` CPUs allocated，约 `<N>` free memory（示例/example: 5/192，约 913 GiB）。 / Slurm `fat`: `<N>` × `<N>`-core node with about `<N>` RAM (example: one 192-core node, about 1.03 TiB); `<N>`/`<N>` CPUs were allocated and about `<N>` memory was free at the snapshot (example: 5/192, about 913 GiB).
- GPU：未发现 `nvidia-smi`，当前分区也未报告 GPU GRES；本阶段按 CPU-only 规划（示例/example: 检查时用户无 Slurm 作业）。 / `nvidia-smi` was not found and current partitions reported no GPU GRES; plan this stage as CPU-only (example: the user had no Slurm jobs at the time).

后续复查 / Later recheck（`<日期 / date>`，示例/example: 2026-08-24）：`free -h` 中的 `<N>` GiB 是 **free swap**，不是可用 RAM（示例/example: 16 GiB）；登录节点 RAM available 约 `<N>` GiB（示例/example: 108 GiB）。项目文件系统可用空间已从约 `<N>` GiB 回升到约 `<N>` GiB（示例/example: 1 → 40），inode 使用率约 `<N>`%（示例/example: 35），但 `df` 仍按整数显示 100% used。QC 表格输出很小，可以运行；大型中间文件仍不应写入项目盘。 / The `<N>` GiB shown by `free -h` is **free swap**, not available RAM (example: 16 GiB); login-node available RAM is about `<N>` GiB (example: 108 GiB). Project-filesystem free space later increased from about `<N>` GiB to about `<N>` GiB (example: 1 → 40) and inode usage is about `<N>`% (example: 35), although integer-rounded `df` still displays 100% used. Small QC tables can be written safely; large intermediates should still not target the project filesystem.
