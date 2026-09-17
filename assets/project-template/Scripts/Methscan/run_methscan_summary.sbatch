#!/usr/bin/env bash
# Validate and summarize all completed threshold branches.
# Usage: sbatch --dependency=afterok:JOB1:JOB2:JOB3 run_methscan_summary.sbatch OUTPUT_DIR

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
[[ -f "$output_dir/methdiff.COMPLETE" ]] || { echo "Missing completed meth-diff stage" >&2; exit 3; }
[[ -f "$output_dir/hypo_heatmaps.COMPLETE" ]] || { echo "Missing completed hypo-DMR heatmap stage" >&2; exit 3; }
for threshold in "${scan_var_thresholds[@]}"; do [[ -f "$output_dir/var_${threshold}.COMPLETE" ]] || { echo "Missing completed branch: $threshold" >&2; exit 3; }; done
args=(--run-dir "$output_dir")
for threshold in "${scan_var_thresholds[@]}"; do args+=(--threshold "$threshold"); done
"$scanpy_python" "$script_dir/05_summarize_run.py" "${args[@]}"
touch "$output_dir/summary.COMPLETE"
printf 'MethSCAn summary completed: %s\n' "$output_dir"
