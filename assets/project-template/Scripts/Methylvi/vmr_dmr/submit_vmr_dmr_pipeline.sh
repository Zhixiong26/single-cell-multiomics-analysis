#!/usr/bin/env bash
# Submit count construction, the VMR+DMR inputs, their MethylVI models, and summary.
# DMR preparation is reused only when its validated COMPLETE marker already exists.
# Stage resources come from config/scheduler.yaml profiles.
set -Eeuo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$script_dir/../../.." && pwd)}
source "$script_dir/00_vmr_dmr_config.sh"
source "$script_dir/../../Common/submit_helpers.sh"
scripts="$SCMO_PROJECT_ROOT/Scripts/Methylvi/vmr_dmr"
root="$SCMO_VMR_DMR_ROOT"
thresholds=("${SCMO_VMR_DMR_THRESHOLDS[@]}")
feature_targets=("${SCMO_VMR_DMR_FEATURE_TARGETS[@]}")
cd "$SCMO_PROJECT_ROOT"

if [[ -f "$SCMO_VMR_DMR_SHARED/prepare.COMPLETE" ]]; then
  prepare_job=reused
  count_job=$(scmo_submit feature_builder "" "$scripts/slurm/build_all_unique_pooled_dmr_counts.sbatch")
else
  [[ ! -e "$root" ]] || { echo "Partial output root exists without prepare.COMPLETE: $root" >&2; exit 2; }
  prepare_job=$(scmo_submit dmr_prepare "" "$scripts/slurm/prepare_all_unique_pooled_dmrs.sbatch")
  count_job=$(scmo_submit feature_builder "afterok:${prepare_job}" \
    "$scripts/slurm/build_all_unique_pooled_dmr_counts.sbatch")
fi

join_jobs=()
train_jobs=()
previous_job=$count_job
for threshold in "${thresholds[@]}"; do
  for features in "${feature_targets[@]}"; do
    join=$(scmo_submit dmr_prepare "afterok:${previous_job}" \
      "$scripts/slurm/combine_vmr_dmr_input.sbatch" "$threshold" "$features")
    train=$(scmo_submit trainer "afterok:${join}" \
      "$scripts/slurm/train_vmr_dmr_model.sbatch" "$threshold" "$features")
    join_jobs+=("$join"); train_jobs+=("$train")
    previous_job=$train
  done
done
summary_job=$(scmo_submit summary "afterok:${previous_job}" "$scripts/slurm/summarize_vmr_dmr_models.sbatch")
printf 'prepare=%s\ncount=%s\njoins=%s\ntrains=%s\nsummary=%s\n' \
  "$prepare_job" "$count_job" "${join_jobs[*]}" "${train_jobs[*]}" "$summary_job"
