#!/usr/bin/env bash
# Submit count construction, six VMR+DMR inputs, six MethylVI models, and summary.
# DMR preparation is reused only when its validated COMPLETE marker already exists.
set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
scripts="$project/Scripts/Methylvi/vmr_dmr"
root="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR"
cd "$project"

if [[ -f "$root/shared/dmr/prepare.COMPLETE" ]]; then
  prepare_job=reused
  count_job=$(sbatch --parsable "$scripts/slurm/build_all_unique_pooled_dmr_counts.sbatch")
else
  [[ ! -e "$root" ]] || { echo "Partial output root exists without prepare.COMPLETE: $root" >&2; exit 2; }
  prepare_job=$(sbatch --parsable "$scripts/slurm/prepare_all_unique_pooled_dmrs.sbatch")
  count_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" \
    "$scripts/slurm/build_all_unique_pooled_dmr_counts.sbatch")
fi

join_jobs=()
train_jobs=()
previous_job=$count_job
for threshold in 0.01 0.02 0.05; do
  for features in 10000 30000; do
    join=$(sbatch --parsable --dependency="afterok:${previous_job}" \
      "$scripts/slurm/combine_vmr_dmr_input.sbatch" "$threshold" "$features")
    train=$(sbatch --parsable --dependency="afterok:${join}" \
      "$scripts/slurm/train_vmr_dmr_model.sbatch" "$threshold" "$features")
    join_jobs+=("$join"); train_jobs+=("$train")
    previous_job=$train
  done
done
summary_job=$(sbatch --parsable --dependency="afterok:${previous_job}" "$scripts/slurm/summarize_vmr_dmr_models.sbatch")
printf 'prepare=%s\ncount=%s\njoins=%s\ntrains=%s\nsummary=%s\n' \
  "$prepare_job" "$count_job" "${join_jobs[*]}" "${train_jobs[*]}" "$summary_job"
