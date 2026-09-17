# Troubleshooting and recovery

- Input signature changed: create a new run ID; do not reuse checkpoints.
- Resource floor cannot fit: wait or choose an allowed larger partition; do not undersize.
- Login path visible but compute path inaccessible: validate from the intended execution context.
- Failure after training: recover only missing plots/summary when signatures match.
- Annotation guard failed: review current markers and clusters; do not reuse old numeric mappings.
- Fractional counts, `mc > cov`, duplicates, or manifest mismatch are hard failures.
- Submitted Slurm jobs are snapshots; later script edits do not alter them.
