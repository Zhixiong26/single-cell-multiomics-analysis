#!/usr/bin/env bash
# Verify all configured model outputs after their dependent training jobs finish.
# Usage: sbatch [resource options] summarize_methylvi_8models.sbatch
# Resources and log paths are injected by the submitter, not declared here.
#SBATCH --job-name=scmo_mvi_summary

set -Eeuo pipefail
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi
source "$script_dir/00_experiment_config.sh"
mvi_root="$SCMO_MVI_EXPERIMENT_ROOT"
summary="$mvi_root/run_summary.tsv"
# The ALLCools route's version column is its 5-kb binning label; the bin size is
# a route parameter, so the label is derived from it rather than written out.
bin_label="$(( ${SCMO_BIN_SIZE:-5000} / 1000 ))kb"
mkdir -p "$mvi_root"
printf 'method\tversion\tfeatures\tstatus\tresults\n' > "$summary"
incomplete=0
check_model() {
  local method=$1 version=$2 features=$3 route=$4
  # Collect what is missing instead of letting `set -e` abort on the first failing test. This
  # wrapper exists to report every configured model, and a model can legitimately lack the
  # supervised UMAP while its clusters are still unannotated. Incompleteness stays fatal below,
  # so workflow.COMPLETE is still never written for a partially produced set.
  local missing=()
  [[ -f "$route/model.COMPLETE" ]] || missing+=(model.COMPLETE)
  [[ -s "$route/results/run_summary.json" ]] || missing+=(results/run_summary.json)
  [[ -s "$route/results/methylvi_embedding.h5ad" ]] || missing+=(results/methylvi_embedding.h5ad)
  [[ -s "$route/results/supervised_umap/supervised_umap_summary.json" ]] || missing+=(results/supervised_umap/supervised_umap_summary.json)
  [[ -s "$route/results/qc/methylation_qc_summary.json" ]] || missing+=(results/qc/methylation_qc_summary.json)
  if (( ${#missing[@]} )); then
    printf '%s\t%s\t%s\tincomplete: %s\t%s\n' "$method" "$version" "$features" "${missing[*]}" "$route/results" >> "$summary"
    incomplete=$(( incomplete + 1 ))
    return 0
  fi
  printf '%s\t%s\t%s\tcomplete\t%s\n' "$method" "$version" "$features" "$route/results" >> "$summary"
}
for features in $SCMO_FEATURE_TARGETS; do
  check_model allcools "$bin_label" "$features" "$mvi_root/allcools/features_${features}"
done
for threshold in $SCMO_VMR_THRESHOLDS; do
  for features in $SCMO_FEATURE_TARGETS; do
    check_model vmr "$threshold" "$features" "$mvi_root/vmr/var_${threshold}/features_${features}"
  done
done
if (( incomplete )); then
  printf 'Incomplete MethylVI models: %d of the configured set; see %s\n' "$incomplete" "$summary" >&2
  exit 1
fi
touch "$mvi_root/workflow.COMPLETE"
printf 'Verified all configured MethylVI models: %s\n' "$mvi_root"
