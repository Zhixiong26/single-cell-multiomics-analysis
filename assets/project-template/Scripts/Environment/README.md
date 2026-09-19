# Environment / 环境

本阶段发现、选择并验证项目所需的软件环境，支持共享环境和分阶段独立环境。

This stage discovers, selects, and validates software environments, supporting both a shared environment and isolated per-stage environments.

## 策略 / Strategy

- 优先复用已有且验证通过的环境。 / Prefer existing verified environments.
- 仅在依赖兼容时共用环境；存在版本冲突时按阶段隔离。 / Share only when dependencies are compatible; isolate stages when versions conflict.
- 未经使用者确认不创建、安装、升级或删除环境。 / Do not create, install, upgrade, or remove environments without user approval.
- 调度脚本使用已记录的绝对可执行路径。 / Scheduler scripts use recorded absolute executable paths.

最终映射维护在 `config/environments.tsv`；探测结果写入 `Results/Environment/`，日志写入 `Scripts/logs/`。

Maintain the final mapping in `config/environments.tsv`; write discovery results to `Results/Environment/` and logs to `Scripts/logs/`.

当前项目采用混合策略：兼容步骤可共享环境，存在依赖冲突的步骤使用独立环境。 / The current project uses a hybrid strategy: compatible stages share an environment, while conflicting stages remain isolated.

环境根目录由使用者/站点配置决定（`$SCMO_ENV_ROOTS`）；发现脚本只读扫描，不修改任何已发现的环境。 / Environment roots come from user/site configuration (`$SCMO_ENV_ROOTS`); the discovery script scans read-only and never modifies a discovered environment.

## 执行入口 / Entry point

```bash
bash Scripts/Environment/01_discover_environments.sh
```

## 环境引导 / Environment bootstrap

发现脚本只发现和登记，不创建环境；缺失的 profile 按 `environment-specs/` 中的版本化 spec 在隔离前缀中创建，不修改共享环境。 / The discovery script only discovers and records; it never creates an environment. Missing profiles are created in isolated prefixes from the versioned specs under `environment-specs/`, without modifying shared environments.

先只读规划，再执行创建： / Plan read-only first, then create:

```bash
PYTHON=/path/to/python-3.9-or-newer
"$PYTHON" tools/bootstrap_environments.py --project "$SCMO_PROJECT_ROOT"
"$PYTHON" tools/bootstrap_environments.py --project "$SCMO_PROJECT_ROOT" --execute
```

默认前缀为项目 `.environments/` 下的 `analysis-core`、`methscan` 和 `methylvi`；profile 的验证命令通过后才登记为可用，创建失败的前缀保留用于诊断而不是静默删除。 / Default prefixes are `analysis-core`, `methscan`, and `methylvi` under the project's `.environments/`; a profile is recorded as ready only after its verification command passes, and a failed creation is retained for diagnosis rather than silently deleted.
