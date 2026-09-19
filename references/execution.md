# Scheduler and execution

`plan_workflow.py` writes a closed DAG under `.workflow/runs/<run-id>/`; writable analysis outputs go to `Results/runs/<run-id>/`. `submit_workflow.py` revalidates signatures, requires matching full validation, and resumes only incomplete tasks.

For Slurm, every task independently queries `sinfo`, `scontrol show nodes`, and `squeue` immediately before `sbatch`. Allocation capacity uses `CPUTot-CPUAlloc` and `RealMemory-AllocMem-headroom`; `CPULoad` and `FreeMem` are diagnostic only. The balanced selector requests the target when it fits, may reduce toward the floor, and fails when no permitted active node meets the floor.

Resource options and dependencies are passed on the `sbatch` command line. Nodes are not pinned unless a profile explicitly requires it. The local backend applies the same resource floors against configured local limits.

Built-in adapters cover every planned task and choose the task's declared stage Python from `environments.tsv`. The task registry is checked during planning and again before production submission. `analysis.task_commands` may override an exact or prefix-wildcard task and may use `{python}`, `{project}`, `{run_id}`, `{run_dir}`, `{result_dir}`, `{task}`, `{task_dir}`, `{threshold}`, and `{features}`.

After a successful override, write `<task_dir>/task_outputs.json` as a JSON object with a non-empty `artifacts` list. Artifact paths should be absolute and every path must exist. The wrapper stores the child command's code as `process_return_code`; missing or invalid evidence makes the wrapper `return_code` nonzero with `failure_stage: output_validation`.

Production stages are entered through `plan_workflow.py` and `submit_workflow.py`. The template also ships standalone `run_*.sbatch` and `submit_*.sh` wrappers for exploring one stage, smoke testing a change, and recovering a single stage of a completed run. They are not a second production path: they produce no `task_outputs.json`, carry no site directive, and take their resources from the same named scheduler profiles the DAG uses, so a wrapper run and a DAG run request identical allocations. A wrapper's success is a stage result, never a run's completion evidence.

The MethylVI route runners (`Methylvi/allcools/run.sh`, `Methylvi/vmr/run.sh`) refuse to start without a managed run ID, because a stage run outside the DAG cannot contribute completion evidence. The packaged wrappers that invoke them set `SCMO_STANDALONE_ACK=1` to acknowledge that explicitly, which is why they work standalone and also print a warning that no completion evidence is recorded. Calling `run.sh` directly still requires either a managed run ID or that acknowledgement; it is never granted implicitly.

Slurm requires a non-empty partition allow-list. Every dependency must exist in the plan; unknown dependencies are fatal. A scheduler `COMPLETED` state is not completion without valid status, signature, marker, and artifacts.

Inspecting a run also refreshes that run's record in the root Report's run log. That write is best effort and never affects the run's verdict: the exit code of `inspect_run.py` reflects the run's state, so an unwritable Report is reported on stderr with the recovery command and leaves the code alone. Use `tools/update_report.py` when the write must fail loudly, or to rebuild the whole log after pruning, migration, or hand edits. See [documentation](documentation.md) for the region format.

Before a run can read its inputs, they have to be where it expects them. Generation links each declared `data.sources` entry into `Data/`, and `tools/link_data.py --project PROJECT --execute` performs the downloads it left pending, verifying each against its declared sha256. Validation refuses any entry that does not resolve from the machine it runs on, so check data reachability in the execution context rather than on the login node. See [configuration](configuration.md) for the block's schema.
