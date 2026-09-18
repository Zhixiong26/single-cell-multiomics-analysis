#!/usr/bin/env bash
# Run pairwise cell-type MethSCAn differential methylation after smooth.
# Usage: sbatch --dependency=afterok:SMOOTH_JOB run_methscan_methdiff.sbatch OUTPUT_DIR

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
[[ -f "$output_dir/smooth.COMPLETE" ]] || { echo "Smooth stage is not complete" >&2; exit 3; }
run_dir="$output_dir/07_methdiff"
"$methscan_python" "$script_dir/07_methdiff_celltype.py" \
  --methscan "$methscan_exe" --data-dir "$output_dir/03_filtered" \
  --cell-metadata "$output_dir/00_scanpy_selected/input_manifest.tsv" \
  --expected-cell-ids "$output_dir/03_filtered/column_header.txt" \
  --chrom-sizes "$SCMO_CHROM_SIZES" \
  --output-dir "$run_dir" --min-cells "$methdiff_min_cells" \
  --threads "$methdiff_threads" --jobs "$methdiff_jobs" \
  --bandwidth "$methdiff_bandwidth" --stepsize "$methdiff_stepsize" \
  --threshold "$methdiff_threshold" --resume
touch "$output_dir/methdiff.COMPLETE"
printf 'Pairwise cell-type meth-diff completed: %s\n' "$run_dir"
