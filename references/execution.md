# Scheduler and execution

`plan_workflow.py` writes a closed DAG under `.workflow/runs/<run-id>/`; writable analysis outputs go to `Results/runs/<run-id>/`. `submit_workflow.py` revalidates signatures, requires matching full validation, and resumes only incomplete tasks.

For Slurm, every task independently queries `sinfo`, `scontrol show nodes`, and `squeue` immediately before `sbatch`. Allocation capacity uses `CPUTot-CPUAlloc` and `RealMemory-AllocMem-headroom`; `CPULoad` and `FreeMem` are diagnostic only. The balanced selector requests the target when it fits, may reduce toward the floor, and fails when no permitted active node meets the floor.

Resource options and dependencies are passed on the `sbatch` command line. Nodes are not pinned unless a profile explicitly requires it.

## Host roles

Where a task may execute is a property of the host, and it is the one policy a config file cannot state: a project that declares `backend: local` while sitting on a cluster login node would still run its tasks on that login node. Every decision point resolves the same role through `_common.host_role()`.

| Role | Decided by | Local execution |
|---|---|---|
| `compute` | `SLURM_JOB_ID`/`SLURM_JOBID` is set: we are inside an allocation | Legitimate — the executor is itself on a compute node |
| `submit` | No allocation, and `scontrol ping` reports a controller `is UP` | Not used: the run is promoted to Slurm |
| `standalone` | No allocation and no controller answers | Legitimate — there is no other host to submit to |

`scontrol ping` is the probe rather than `sinfo` or `squeue`. Slurm client commands on PATH prove nothing, `squeue -u $USER` prints nothing and exits 0 on a laptop and on a login node alike, and a controller that answers but is `DOWN` refuses every submission — so the reply itself has to report UP, not merely arrive.

A `local` project on a submit host is **promoted to Slurm**, not refused and not executed in place. The promotion is recorded: `submissions.json` carries `backend_declared` beside `backend_effective`, each `resource_snapshots/*.json` carries `host_role`, and every `task_status.json` records the role and `slurm_job_id` the task saw. A promotion needs a partition allow-list; when `scheduler.partitions` is empty it refuses and names the partitions `sinfo` reports.

**A planned task is never executed on a submit host.** The guard is at execution rather than at planning: `Scripts/Common/run_task.py` refuses before it reads a plan or creates a task directory, and `Scripts/Common/task_adapter.py` refuses before it runs a stage, because a hand call to either is a separate entry point. `Methylvi/*/run.sh` carry the same refusal in shell. A refusal is a precondition failure, not a task result: it exits 2, writes no `task_status.json`, and cannot be mistaken for a failed attempt. The legitimate work of a login node is planning, validation, submission, inspection, and reporting — including the `sinfo`/`squeue` queries the resource check performs.

The one exception is explicit and run-scoped. `submit_workflow.py --allow-login-execution` names a single run, is passed to the task as `SCMO_LOGIN_EXECUTION_ACK=<run-id>`, and is honoured only when that value equals the run being executed, so it cannot carry over to the next run. It is never inherited by a submitted job or written into `os.environ` — an inherited bypass would disarm the guard in the very child process it exists to constrain. Every task it covers is recorded as `login_execution: true`. Prefer the alternative it stands in for: run the orchestrator inside an allocation (`salloc`, or `srun --pty bash`), where the role is `compute` and the local executor is on a compute node.

Resource inspection for the `local` backend is against the local machine: the recommendation is the profile target clamped to `scheduler.local.max_threads`/`max_memory` and floored against the profile, and no cluster is queried, because none is involved. A snapshot whose `backend` is `local` therefore reports `nodes: []` by design, and says so in its `note`.

Built-in adapters cover every planned task and choose the task's declared stage Python from `environments.tsv`. The task registry is checked during planning and again before production submission. `analysis.task_commands` may override an exact or prefix-wildcard task and may use `{python}`, `{project}`, `{run_id}`, `{run_dir}`, `{result_dir}`, `{task}`, `{task_dir}`, `{threshold}`, and `{features}`.

After a successful override, write `<task_dir>/task_outputs.json` as a JSON object with a non-empty `artifacts` list. Artifact paths should be absolute and every path must exist. The wrapper stores the child command's code as `process_return_code`; missing or invalid evidence makes the wrapper `return_code` nonzero with `failure_stage: output_validation`.

Production stages are entered through `plan_workflow.py` and `submit_workflow.py`. The template also ships standalone `run_*.sbatch` and `submit_*.sh` wrappers for exploring one stage, smoke testing a change, and recovering a single stage of a completed run. They are not a second production path: they produce no `task_outputs.json`, carry no site directive, and take their resources from the same named scheduler profiles the DAG uses, so a wrapper run and a DAG run request identical allocations. A wrapper's success is a stage result, never a run's completion evidence.

The MethylVI route runners (`Methylvi/allcools/run.sh`, `Methylvi/vmr/run.sh`) refuse to start without a managed run ID, because a stage run outside the DAG cannot contribute completion evidence. The packaged wrappers that invoke them set `SCMO_STANDALONE_ACK=1` to acknowledge that explicitly, which is why they work standalone and also print a warning that no completion evidence is recorded. Calling `run.sh` directly still requires either a managed run ID or that acknowledgement; it is never granted implicitly. That acknowledgement says a run may bypass the DAG and says nothing about which host may execute it, so on a submit host a direct call is refused as well.

The packaged `run_*.sbatch` wrappers are sbatch entry points: `submit_*.sh` submits them, and the DAG reaches them from inside an allocation. They are `bash` scripts with no execution guard of their own, so `bash Scripts/Methscan/run_methscan_common.sbatch OUT` on a login node would run that stage in the caller's process. Submit them; do not source or bash them by hand outside an allocation.

Slurm requires a non-empty partition allow-list. Every dependency must exist in the plan; unknown dependencies are fatal. A scheduler `COMPLETED` state is not completion without valid status, signature, marker, and artifacts.

Inspecting a run also refreshes that run's record in the root Report's run log. That write is best effort and never affects the run's verdict: the exit code of `inspect_run.py` reflects the run's state, so an unwritable Report is reported on stderr with the recovery command and leaves the code alone. Use `tools/update_report.py` when the write must fail loudly, or to rebuild the whole log after pruning, migration, or hand edits. See [documentation](documentation.md) for the region format.

Before a run can read its inputs, they have to be where it expects them. Generation links each declared `data.sources` entry into `Data/`, and `tools/link_data.py --project PROJECT --execute` performs the downloads it left pending, verifying each against its declared sha256. Validation refuses any entry that does not resolve from the machine it runs on, so check data reachability in the execution context rather than on the login node. See [configuration](configuration.md) for the block's schema.

Provisioning and integrity work is not a planned task and is not covered by the execution guard: `validate_project.py --mode full` reads every ALLC file end to end, `bootstrap_environments.py --execute` solves and creates environments, and `link_data.py --execute` transfers downloads. Each runs where it is invoked, because each has to — an environment belongs in the prefix it creates, and a data check has to see the filesystem the run will see. On a site that forbids such work on a login node, run those three inside an allocation. Their results are not a run's completion evidence either way.
