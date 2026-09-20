---
name: single-cell-multiomics-analysis
description: Generate, validate, schedule, monitor, and document portable single-cell RNA and DNA-methylation projects using Scanpy/Harmony, MethSCAn, ALLCools, and MethylVI. Use for configuration-driven project creation or execution from 10x RNA matrices and per-cell indexed ALLC data; do not use for raw FASTQ/BAM processing.
---

# Single-cell multiomics analysis

Create auditable projects from user-supplied data without assuming sample names, paths, cell counts, genome builds, environments, or cluster annotations. Use the bundled generator and templates instead of copying paths from a previous project.

## Operating contract

1. Read [intake and configuration](references/configuration.md), then create an intake YAML with explicit sample metadata, inputs, genome resources, any known environments, and scheduler site policy. Declare where the user's data already lives under `data.sources`; do not point `config/` at a path outside the project.
2. Generate a new project with `scripts/init_project.py`. Never generate into a non-empty directory. Generation links each declared source into the project's own `Data/`; a `url` source is fetched by `tools/link_data.py --execute`, which verifies its sha256.
3. Read [environment provisioning](references/environments.md). Discover compatible existing environments first. If a required environment is unavailable, use `tools/bootstrap_environments.py --execute` to create a new isolated project environment, verify it, and update `config/environments.tsv`.
4. Run quick validation before planning and immediately before submission. A matching full-validation record is required for first production use or changed inputs/schema/template. Missing inputs, indexes, required references, declared-checksum matches, route-used environments, or route prerequisites are blockers; an explicit route subset must not require unrelated stages.
5. Read the relevant workflow reference before adapting task commands: [Scanpy](references/scanpy.md), [MethSCAn](references/methscan.md), or [ALLCools/MethylVI](references/methylvi.md). A Scanpy pass is not finished when the clusters are still `Unassigned`: read that run's `candidate_annotation_audit.tsv`, establish each cluster's cell type from its own ranked markers looked up as you go, record the mapping with `tools/record_annotation_review.py`, re-run as a baseline, and then ask the user to confirm the labels you are sure of and the ones you are not. See [Scanpy](references/scanpy.md) for the loop.
6. Build a distinct run DAG with `scripts/plan_workflow.py`. Results belong under `Results/runs/<run-id>/`; reuse a run ID only to resume its matching signed plan.
7. For execution, read [scheduler and execution](references/execution.md). Submit only when the user requested execution. Every task must pass a fresh resource inspection immediately before its `sbatch`; a prior snapshot is not reusable. Never execute a planned task on a login node: inspect the compute nodes' resources and submit, and let the login node do only planning, validation, submission, inspection, and reporting.
8. Judge completion from validated outputs, task completion markers, and `run_summary.json`, not from submission success or Slurm state alone.
9. Keep the root README as the stable contract: `init_project.py` writes it once from the intake and the resolved configuration, and nothing ever rewrites it with runtime state. Every `inspect_run.py` refreshes that run's concise record in the root Report's run log *and* in each touched stage Report's own log, so the record a reader of `Scripts/Scanpy/Report.md` sees is the same evidence filtered to Scanpy's tasks. `update_report.py` rebuilds the root log and every stage log from `.workflow/runs/*/run_summary.json` after pruning, migration, or hand edits. The stage READMEs and the prose around the stage logs are generated-once text; `init_project.py` binds them to the project at generation and no run rewrites them.

## Invariants

- Treat RNA, ALLC, annotation, genome, and pre-existing shared environments as read-only inputs.
- Every analysis input must be reachable from the project as `Data/<name>`. Link or download it there rather than configuring an absolute path elsewhere; a project that reads outside its own directory is not reproducible from it, and the sample manifest is where the deviation becomes invisible.
- Never upgrade or mutate a discovered shared environment. Create missing dependencies under the generated project's `.environments/` directory or an explicitly configured new prefix.
- Do not infer condition, batch, donor, tissue, or cell identity from directory names.
- Require `.tbi` beside every `.allc.tsv.gz` used as indexed ALLC input.
- Derive ALLC ownership from `samples.tsv`, chromosomes from `chrom_sizes`, and methylation context from configuration; never infer them from a previous project.
- Preserve integer MethylVI `mc`/`cov`, require `0 <= mc <= cov`, unique cells/features, and matching input signatures before checkpoint reuse.
- A new Scanpy cluster set is unreviewed. Never transplant numeric cluster mappings without a matching analysis signature and exact cluster set; `tools/record_annotation_review.py` reviews the run's own clusters and writes both into the profile it records.
- Annotate every cluster in the pass that produced it, from that run's own marker evidence. Do not leave a cluster at `Unassigned` as a resting state, and do not use a bundled marker panel or another project's crosswalk as the source of a label; a label you cannot defend from the run's markers stays explicitly uncertain and is reported as such. The pass is done when the delivered annotation dotplot has one row per cell type with the samples pooled into it: rows are cell types, never clusters and never samples, and a per-cluster fallback dotplot still standing alone means the labels were not settled.
- Never use `NA`, `Unassigned`, `requires_review`, or an unapproved annotation table for cell-type DMR.
- Keep Top-N hypo-DMR heatmap selection separate from the all-unique pooled-DMR MethylVI route.
- Never execute a planned task on a Slurm submit host. Resource inspection and submission belong there; execution belongs on a compute node or, on a machine with no scheduler at all, on that machine. `run_task.py` and `task_adapter.py` refuse at execution, and `Methylvi/*/run.sh` refuse in shell; a refusal is a precondition failure that writes no task evidence, not a failed task. If work must run on the machine you are on, use an allocation (`salloc`, `srun --pty bash`) rather than an override.
- Execute production tasks only through the planned DAG. Overrides must emit `task_outputs.json` with a non-empty list of existing artifacts. The packaged `run_*.sbatch` and `submit_*.sh` wrappers are for exploration, smoke testing, and recovery; they resolve resources from scheduler profiles and are never the source of a run's completion evidence. The `.sbatch` wrappers are sbatch entry points — submit them rather than running them by hand outside an allocation.
- Report Slurm-reserved CPU/memory separately from observed load/`FreeMem`. Never treat reserved memory as proof of actual use.
- Do not reduce below a task resource floor merely to fit a busy node. Do not pin a node unless data visibility or hardware requirements demand it.
- Use a small isolated test when the template, environment, schema, or high-risk route changed; ordinary repeated runs still require full read-only preflight.

## Packaged project code

`assets/project-template/` contains portable copies of the maintained Environment, Scanpy, MethSCAn, ALLCools, MethylVI, VMR, and VMR+DMR code. Generated projects receive their own runtime tools and remain runnable without this skill installed. Read [template maintenance](references/template-maintenance.md) before updating packaged algorithms.

## Conditional references

- Read [troubleshooting](references/troubleshooting.md) for failed jobs, inaccessible storage, signature mismatch, or recovery.
- Read [documentation](references/documentation.md) when creating or updating README, Report, manifests, and completion evidence.
- Read [Git collaboration](references/git.md) only when the user asks to initialize, commit, connect, or push a repository.
