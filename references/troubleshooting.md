# Troubleshooting and recovery

- Input signature changed: create a new run ID; do not reuse checkpoints. Reuse the same run ID only to resume its unchanged signed plan.
- Missing full-validation evidence: run `validate_project.py --mode full`; production submit also creates and caches matching evidence when absent.
- Missing task artifact with a stale marker: resubmit the same run ID; the task is treated as incomplete and rerun.
- Resource floor cannot fit: wait or choose an allowed larger partition; do not undersize.
- Login path visible but compute path inaccessible: validate from the intended execution context.
- Failure after training: recover only missing plots/summary when signatures match.
- Annotation guard failed: review current markers and clusters; do not reuse old numeric mappings.
- DMR route disabled: set `annotation.review_status: approved` only after reviewing labels and remove placeholder cell types.
- Fractional counts, `mc > cov`, duplicates, or manifest mismatch are hard failures.
- Submitted Slurm jobs are snapshots; later script edits do not alter them.
- Missing environment stage: run `tools/bootstrap_environments.py --project PROJECT` to inspect the plan, then add `--execute`. Do not repair it by upgrading a shared environment.
- Override command exited zero but task failed: inspect `process_return_code`, `failure_stage`, and `error` in `task_status.json`; create `<task_dir>/task_outputs.json` with a non-empty list of existing artifacts.
- Failed Conda creation: inspect `.workflow/environment-bootstrap/` and the incomplete isolated prefix. The bootstrapper intentionally does not delete it; choose a new empty prefix after correcting channel, network, or storage problems.
