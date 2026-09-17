#!/usr/bin/env bash
# Pool configured samples by cell type and run all pairwise MethSCAn comparisons.
# Usage: sbatch run_methscan_pooled_methdiff.sbatch OUTPUT_DIR

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
[[ -f "$output_dir/smooth.COMPLETE" ]] || { echo "Smooth stage is not complete" >&2; exit 3; }
pooled_root="$output_dir/pooled"
mkdir -p "$pooled_root"
"$methscan_python" "$script_dir/07_methdiff_celltype.py" \
  --methscan "$methscan_exe" --data-dir "$output_dir/03_filtered" \
  --cell-metadata "$output_dir/00_scanpy_selected/input_manifest.tsv" \
  --expected-cell-ids "$output_dir/03_filtered/column_header.txt" \
  --output-dir "$pooled_root/07_methdiff" --min-cells "$methdiff_min_cells" \
  --threads "$methdiff_threads" --jobs "$methdiff_jobs" \
  --bandwidth "$methdiff_bandwidth" --stepsize "$methdiff_stepsize" \
  --threshold "$methdiff_threshold" --pool-samples-as "$pooled_sample_label" --resume
touch "$pooled_root/methdiff.COMPLETE"
printf 'Pooled pairwise cell-type meth-diff completed: %s\n' "$pooled_root/07_methdiff"
