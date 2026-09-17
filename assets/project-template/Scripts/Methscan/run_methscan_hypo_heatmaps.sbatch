#!/usr/bin/env bash
# Resolve optional raw-p fallback, build Top200 hypo-DMR matrices, and plot.
# Usage: sbatch --dependency=afterok:METHDIFF_JOB run_methscan_hypo_heatmaps.sbatch OUTPUT_DIR

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
[[ -f "$output_dir/methdiff.COMPLETE" ]] || { echo "Meth-diff stage is not complete" >&2; exit 3; }
"$scanpy_python" "$script_dir/08_hypo_dmr_heatmaps.py" \
  --run-dir "$output_dir" --methscan "$methscan_exe" --methscan-python "$methscan_python" \
  --raw-p "$hypo_raw_p" --min-abs-diff "$hypo_min_abs_diff" \
  --top-per-cell-type "$hypo_top_per_cell_type" --min-cells "$methdiff_min_cells" \
  --threads "$methdiff_threads" --bandwidth "$methdiff_bandwidth" \
  --stepsize "$methdiff_stepsize" --threshold "$methdiff_threshold" \
  --matrix-workers "${SLURM_CPUS_PER_TASK:-$hypo_matrix_workers}" \
  --zscore-min-observed-cells "$hypo_zscore_min_observed_cells" \
  --zscore-clip "$hypo_zscore_clip" \
  --compressed-colorbar-limit "$hypo_compressed_colorbar_limit" \
  --dpi "$hypo_heatmap_dpi"
touch "$output_dir/hypo_heatmaps.COMPLETE"
printf 'Top200 hypo-DMR matrices and heatmaps completed: %s\n' "$output_dir"
