---
name: single-cell-multiomics-analysis
description: Generate, validate, schedule, monitor, and document portable single-cell RNA and DNA-methylation projects using Scanpy/Harmony, MethSCAn, ALLCools, and MethylVI. Use for configuration-driven project creation or execution from 10x RNA matrices and per-cell indexed ALLC data; do not use for raw FASTQ/BAM processing.
---

# Single-cell multiomics analysis

Create auditable projects from user-supplied data without assuming sample names, paths, cell counts, genome builds, environments, or cluster annotations. Use the bundled generator and templates instead of copying paths from a previous project.

## Operating contract

1. Read [intake and configuration](references/configuration.md), then create an intake YAML with explicit sample metadata, inputs, genome resources, any known environments, and scheduler site policy.
2. Generate a new project with `scripts/init_project.py`. Never generate into a non-empty directory.
3. Read [environment provisioning](references/environments.md). Discover compatible existing environments first. If a required environment is unavailable, use `tools/bootstrap_environments.py --execute` to create a new isolated project environment, verify it, and update `config/environments.tsv`.
4. Run quick validation before planning and immediately before submission. A matching full-validation record is required for first production use or changed inputs/schema/template. Missing inputs, indexes, required references, declared-checksum matches, route-used environments, or route prerequisites are blockers; an explicit route subset must not require unrelated stages.
5. Read the relevant workflow reference before adapting task commands: [Scanpy](references/scanpy.md), [MethSCAn](references/methscan.md), or [ALLCools/MethylVI](references/methylvi.md).
6. Build a distinct run DAG with `scripts/plan_workflow.py`. Results belong under `Results/runs/<run-id>/`; reuse a run ID only to resume its matching signed plan.
7. For execution, read [scheduler and execution](references/execution.md). Submit only when the user requested execution. Every task must pass a fresh resource inspection immediately before its `sbatch`; a prior snapshot is not reusable.
8. Judge completion from validated outputs, task completion markers, and `run_summary.json`, not from submission success or Slurm state alone.
9. Keep README as the stable contract. Update bilingual Report evidence with `inspect_run.py` and `update_report.py` after material progress, failure, or recovery.

## Invariants

- Treat RNA, ALLC, annotation, genome, and pre-existing shared environments as read-only inputs.
- Never upgrade or mutate a discovered shared environment. Create missing dependencies under the generated project's `.environments/` directory or an explicitly configured new prefix.
- Do not infer condition, batch, donor, tissue, or cell identity from directory names.
- Require `.tbi` beside every `.allc.tsv.gz` used as indexed ALLC input.
- Derive ALLC ownership from `samples.tsv`, chromosomes from `chrom_sizes`, and methylation context from configuration; never infer them from a previous project.
- Preserve integer MethylVI `mc`/`cov`, require `0 <= mc <= cov`, unique cells/features, and matching input signatures before checkpoint reuse.
- A new Scanpy cluster set is unreviewed. Never transplant numeric cluster mappings without a matching analysis signature and exact cluster set.
- Never use `NA`, `Unassigned`, `requires_review`, or an unapproved annotation table for cell-type DMR.
- Keep Top-N hypo-DMR heatmap selection separate from the all-unique pooled-DMR MethylVI route.
- Execute production tasks only through the planned DAG. Overrides must emit `task_outputs.json` with a non-empty list of existing artifacts; do not run stage-local Slurm wrappers directly.
- Report Slurm-reserved CPU/memory separately from observed load/`FreeMem`. Never treat reserved memory as proof of actual use.
- Do not reduce below a task resource floor merely to fit a busy node. Do not pin a node unless data visibility or hardware requirements demand it.
- Use a small isolated test when the template, environment, schema, or high-risk route changed; ordinary repeated runs still require full read-only preflight.

## Packaged project code

`assets/project-template/` contains portable copies of the maintained Environment, Scanpy, MethSCAn, ALLCools, MethylVI, VMR, and VMR+DMR code. Generated projects receive their own runtime tools and remain runnable without this skill installed. Read [template maintenance](references/template-maintenance.md) before updating packaged algorithms.

## Conditional references

- Read [troubleshooting](references/troubleshooting.md) for failed jobs, inaccessible storage, signature mismatch, or recovery.
- Read [documentation](references/documentation.md) when creating or updating README, Report, manifests, and completion evidence.
- Read [Git collaboration](references/git.md) only when the user asks to initialize, commit, connect, or push a repository.
