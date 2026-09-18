# Environment / 环境

Run `01_discover_environments.sh` read-only. Reuse compatible shared environments without changing them. If a required profile is missing, run `tools/bootstrap_environments.py --project PROJECT --execute` to create and verify an isolated environment under `PROJECT/.environments/` and update `config/environments.tsv`.

先只读发现并复用兼容环境，不修改共享环境；缺失环境通过 `bootstrap_environments.py` 在项目 `.environments/` 下新建、验证并登记。
