#!/usr/bin/env bash
# One-click submission of two ALLCools and six VMR MethylVI models.
# Stage resources come from config/scheduler.yaml profiles. The managed run DAG
# (workflow.py plan/submit) is the production path; this wrapper is the
# equivalent standalone entry point.
set -Eeuo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
source "$script_dir/00_experiment_config.sh"
source "$script_dir/../Common/submit_helpers.sh"
scripts="$SCMO_PROJECT_ROOT/Scripts/Methylvi"
[[ ! -e "$SCMO_MVI_EXPERIMENT_ROOT" ]] || { echo "Experiment root already exists: $SCMO_MVI_EXPERIMENT_ROOT" >&2; exit 2; }

read -r -a thresholds <<< "${SCMO_VMR_THRESHOLDS:?set SCMO_VMR_THRESHOLDS}"
read -r -a feature_targets <<< "${SCMO_FEATURE_TARGETS:?set SCMO_FEATURE_TARGETS}"

allc_build=$(scmo_submit feature_builder "" "$scripts/allcools/slurm/build_allcools_feature_inputs.sbatch")
vmr_builds=()
train_jobs=()
for threshold in "${thresholds[@]}"; do
  build=$(scmo_submit feature_builder "" "$scripts/vmr/slurm/build_vmr_feature_inputs.sbatch" "$threshold")
  vmr_builds+=("$build")
  for features in "${feature_targets[@]}"; do
    train_jobs+=("$(scmo_submit trainer "afterok:${build}" "$scripts/slurm/train_methylvi_model.sbatch" vmr "$threshold" "$features")")
  done
done
for features in "${feature_targets[@]}"; do
  train_jobs+=("$(scmo_submit trainer "afterok:${allc_build}" "$scripts/slurm/train_methylvi_model.sbatch" allcools none "$features")")
done
dependency=$(IFS=:; echo "afterok:${train_jobs[*]}")
summary_job=$(scmo_submit summary "$dependency" "$scripts/slurm/summarize_methylvi_8models.sbatch")
printf 'allcools_build=%s\nvmr_builds=%s\ntraining_jobs=%s\nsummary=%s\n' "$allc_build" "${vmr_builds[*]}" "${train_jobs[*]}" "$summary_job"
