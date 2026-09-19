#!/usr/bin/env bash
# One-click submit: common stages -> three parallel threshold branches -> summary.
# Stage resources come from config/scheduler.yaml profiles; no site value is
# hardcoded here. The managed run DAG (workflow.py plan/submit) is the production
# path; this wrapper is the equivalent standalone entry point.
# Usage: submit_methscan_pipeline.sh OUTPUT_DIR
set -Eeuo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
source "$script_dir/00_methscan_config.sh"
source "$script_dir/../Common/submit_helpers.sh"
if (( $# != 1 )); then echo "Usage: $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1

common_job=$(scmo_submit io_builder "" "$script_dir/run_methscan_common.sbatch" "$output_dir")
prepare_job=$(scmo_submit serial "afterok:${common_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" prepare)
filter_job=$(scmo_submit serial "afterok:${prepare_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" filter)
smooth_job=$(scmo_submit serial "afterok:${filter_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" smooth)
methdiff_job=$(scmo_submit dmr "afterok:${smooth_job}" "$script_dir/run_methscan_methdiff.sbatch" "$output_dir")
hypo_heatmap_job=$(scmo_submit plot "afterok:${methdiff_job}" "$script_dir/run_methscan_hypo_heatmaps.sbatch" "$output_dir")
branch_jobs=()
for threshold in "${scan_var_thresholds[@]}"; do
  branch_jobs+=("$(scmo_submit methscan_branch "afterok:${smooth_job}" "$script_dir/run_methscan_branch.sbatch" "$output_dir" "$threshold")")
done
dependency=$(IFS=:; echo "afterok:${hypo_heatmap_job}:${branch_jobs[*]}")
summary_job=$(scmo_submit summary "$dependency" "$script_dir/run_methscan_summary.sbatch" "$output_dir")
printf 'common=%s\nprepare=%s\nfilter=%s\nsmooth=%s\nmethdiff=%s\nhypo_heatmaps=%s\nbranches=%s\nsummary=%s\n' "$common_job" "$prepare_job" "$filter_job" "$smooth_job" "$methdiff_job" "$hypo_heatmap_job" "${branch_jobs[*]}" "$summary_job"
