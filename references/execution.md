# Scheduler and execution

`plan_workflow.py` writes a closed DAG under `.workflow/runs/<run-id>/`; writable analysis outputs go to `Results/runs/<run-id>/`. `submit_workflow.py` revalidates signatures, requires matching full validation, and resumes only incomplete tasks.

For Slurm, every task independently queries `sinfo`, `scontrol show nodes`, and `squeue` immediately before `sbatch`. Allocation capacity uses `CPUTot-CPUAlloc` and `RealMemory-AllocMem-headroom`; `CPULoad` and `FreeMem` are diagnostic only. The balanced selector requests the target when it fits, may reduce toward the floor, and fails when no permitted active node meets the floor.

Resource options and dependencies are passed on the `sbatch` command line. Nodes are not pinned unless a profile explicitly requires it. The local backend applies the same resource floors against configured local limits.

Built-in adapters cover every planned task and choose the task's declared stage Python from `environments.tsv`. The task registry is checked during planning and again before production submission. `analysis.task_commands` may override an exact or prefix-wildcard task and may use `{python}`, `{project}`, `{run_id}`, `{run_dir}`, `{result_dir}`, `{task}`, `{task_dir}`, `{threshold}`, and `{features}`.

After a successful override, write `<task_dir>/task_outputs.json` as a JSON object with a non-empty `artifacts` list. Artifact paths should be absolute and every path must exist. The wrapper stores the child command's code as `process_return_code`; missing or invalid evidence makes the wrapper `return_code` nonzero with `failure_stage: output_validation`. Production stages must be entered through `plan_workflow.py` and `submit_workflow.py`, not by directly invoking stage-local Slurm files.

Slurm requires a non-empty partition allow-list. Every dependency must exist in the plan; unknown dependencies are fatal. A scheduler `COMPLETED` state is not completion without valid status, signature, marker, and artifacts.
