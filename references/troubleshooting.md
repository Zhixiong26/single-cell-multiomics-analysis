# Troubleshooting and recovery

- Input signature changed: create a new run ID; do not reuse checkpoints.
- Resource floor cannot fit: wait or choose an allowed larger partition; do not undersize.
- Login path visible but compute path inaccessible: validate from the intended execution context.
- Failure after training: recover only missing plots/summary when signatures match.
- Annotation guard failed: review current markers and clusters; do not reuse old numeric mappings.
- Fractional counts, `mc > cov`, duplicates, or manifest mismatch are hard failures.
- Submitted Slurm jobs are snapshots; later script edits do not alter them.
- Missing environment stage: run `tools/bootstrap_environments.py --project PROJECT` to inspect the plan, then add `--execute`. Do not repair it by upgrading a shared environment.
- Failed Conda creation: inspect `.workflow/environment-bootstrap/` and the incomplete isolated prefix. The bootstrapper intentionally does not delete it; choose a new empty prefix after correcting channel, network, or storage problems.
