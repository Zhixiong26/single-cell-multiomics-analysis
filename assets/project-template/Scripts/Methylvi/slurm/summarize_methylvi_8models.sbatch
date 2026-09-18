#!/usr/bin/env bash
# Verify all eight model outputs after their dependent training jobs finish.

set -Eeuo pipefail
project=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
source "$project/Scripts/Methylvi/00_experiment_config.sh"
summary="$MVI_EXPERIMENT_ROOT/run_summary.tsv"
mkdir -p "$MVI_EXPERIMENT_ROOT"
printf 'method\tversion\tfeatures\tstatus\tresults\n' > "$summary"
check_model() {
  local method=$1 version=$2 features=$3 route=$4
  [[ -f "$route/model.COMPLETE" ]]
  [[ -s "$route/results/run_summary.json" ]]
  [[ -s "$route/results/methylvi_embedding.h5ad" ]]
  [[ -s "$route/results/supervised_umap/supervised_umap_summary.json" ]]
  [[ -s "$route/results/qc/methylation_qc_summary.json" ]]
  printf '%s\t%s\t%s\tcomplete\t%s\n' "$method" "$version" "$features" "$route/results" >> "$summary"
}
for features in $MVI_FEATURE_TARGETS; do
  check_model allcools 5kb "$features" "$MVI_EXPERIMENT_ROOT/allcools/features_${features}"
done
for threshold in $MVI_VMR_THRESHOLDS; do
  for features in $MVI_FEATURE_TARGETS; do
    check_model vmr "$threshold" "$features" "$MVI_EXPERIMENT_ROOT/vmr/var_${threshold}/features_${features}"
  done
done
touch "$MVI_EXPERIMENT_ROOT/workflow.COMPLETE"
printf 'Verified configured MethylVI models: %s\n' "$MVI_EXPERIMENT_ROOT"
