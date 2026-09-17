# Scheduler and execution

`plan_workflow.py` writes `.workflow/runs/<run-id>/plan.json` and preflight evidence without submitting. `submit_workflow.py` revalidates the input signature and refuses changed inputs.

For Slurm, every task independently queries `sinfo`, `scontrol show nodes`, and `squeue` immediately before `sbatch`. Allocation capacity uses `CPUTot-CPUAlloc` and `RealMemory-AllocMem-headroom`; `CPULoad` and `FreeMem` are diagnostic only. The balanced selector requests the target when it fits, may reduce toward the floor, and fails when no permitted active node meets the floor.

Resource options and dependencies are passed on the `sbatch` command line. Nodes are not pinned unless a profile explicitly requires it. The local backend applies the same resource floors against configured local limits.

Before execution, ensure each planned task has an exact or prefix-wildcard `analysis.task_commands` entry. Commands may use `{python}`, `{project}`, `{run_id}`, `{run_dir}`, `{task}`, `{task_dir}`, `{threshold}`, and `{features}` placeholders. Use `{python}` instead of an unqualified `python3` to preserve the verified environment.
