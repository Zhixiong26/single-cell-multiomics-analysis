#!/usr/bin/env bash
# One-click submission of two ALLCools and six VMR MethylVI models.
set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
source "$project/Scripts/Methylvi/00_experiment_config.sh"
scripts="$project/Scripts/Methylvi"
[[ ! -e "$MVI_EXPERIMENT_ROOT" ]] || { echo "Experiment root already exists: $MVI_EXPERIMENT_ROOT" >&2; exit 2; }
allc_build=$(sbatch --parsable "$scripts/allcools/slurm/build_allcools_feature_inputs.sbatch")
vmr_builds=()
train_jobs=()
for threshold in 0.01 0.02 0.05; do
  build=$(sbatch --parsable "$scripts/vmr/slurm/build_vmr_feature_inputs.sbatch" "$threshold")
  vmr_builds+=("$build")
  for features in 10000 30000; do
    train_jobs+=("$(sbatch --parsable --dependency="afterok:${build}" "$scripts/slurm/train_methylvi_model.sbatch" vmr "$threshold" "$features")")
  done
done
for features in 10000 30000; do
  train_jobs+=("$(sbatch --parsable --dependency="afterok:${allc_build}" "$scripts/slurm/train_methylvi_model.sbatch" allcools none "$features")")
done
dependency=$(IFS=:; echo "afterok:${train_jobs[*]}")
summary_job=$(sbatch --parsable --dependency="$dependency" "$scripts/slurm/summarize_methylvi_8models.sbatch")
printf 'allcools_build=%s\nvmr_builds=%s\ntraining_jobs=%s\nsummary=%s\n' "$allc_build" "${vmr_builds[*]}" "${train_jobs[*]}" "$summary_job"
