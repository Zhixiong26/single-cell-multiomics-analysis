#!/usr/bin/env bash
# Run one independent VMR threshold branch after the common stages.
# Usage: sbatch run_methscan_branch.sbatch OUTPUT_DIR THRESHOLD

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 2 )); then echo "Usage: sbatch $0 OUTPUT_DIR THRESHOLD" >&2; exit 2; fi
output_dir=$1
threshold=$2
[[ "$(readlink -m "$output_dir")" == "$(readlink -m "${SCMO_METHSCAN_RUN_DIR:?set SCMO_METHSCAN_RUN_DIR}")" ]] || { echo "Unexpected MethSCAn output root" >&2; exit 3; }
[[ " ${scan_var_thresholds[*]} " == *" $threshold "* ]] || { echo "Unsupported threshold: $threshold" >&2; exit 4; }
[[ -f "$output_dir/smooth.COMPLETE" ]] || { echo "Smooth stage is not complete: $output_dir" >&2; exit 5; }
[[ -r "$output_dir/03_filtered/column_header.txt" ]] || { echo "Filtered input is missing" >&2; exit 6; }
[[ -x "$methscan_exe" && -x "$scanpy_python" ]] || { echo "Configured executable missing" >&2; exit 7; }
label="var_${threshold}"
scan_dir="$output_dir/04_scan/$label"
matrix_dir="$output_dir/05_matrix/$label"
scanpy_dir="$output_dir/06_scanpy/$label"
status_file="$output_dir/${label}.status.tsv"
mkdir -p "$scan_dir" "$matrix_dir" "$scanpy_dir"
printf 'stage\tstatus\ttimestamp\n' > "$status_file"
record_stage() { printf '%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" >> "$status_file"; }
current_stage=initialization
on_error() { local code=$?; trap - ERR; record_stage "$current_stage" "failed(exit=$code)"; exit "$code"; }
trap on_error ERR
run_stage() { current_stage=$1; shift; record_stage "$current_stage" running; "$methscan_python" "$script_dir/06_run_with_resources.py" --output "$output_dir/${label}.${current_stage}.resources.json" -- "$@"; record_stage "$current_stage" complete; }
run_stage scan "$methscan_exe" scan --bandwidth "$scan_bandwidth" --stepsize "$scan_stepsize" --var-threshold "$threshold" --min-cells "$scan_min_cells" --threads "${SLURM_CPUS_PER_TASK:-$threads}" "$output_dir/03_filtered" "$scan_dir/VMRs.bed"
run_stage matrix "$methscan_exe" matrix --threads "${SLURM_CPUS_PER_TASK:-$threads}" "$scan_dir/VMRs.bed" "$output_dir/03_filtered" "$matrix_dir"
run_stage scanpy "$scanpy_python" "$script_dir/04_vmr_scanpy.py" --matrix "$matrix_dir/mean_shrunken_residuals.csv.gz" --expected-cell-ids "$output_dir/03_filtered/column_header.txt" --cell-metadata "$output_dir/00_scanpy_selected/input_manifest.tsv" --output-dir "$scanpy_dir" --min-region-cell-fraction "$vmr_min_cell_fraction" --min-cell-regions "$cell_min_regions" --n-pcs "$pca_components" --neighbors "$neighbors" --leiden-resolution "$leiden_resolution" --seed "$random_seed"
touch "$output_dir/${label}.COMPLETE"
record_stage branch workflow_complete
printf 'Threshold branch completed: %s (%s)\n' "$output_dir" "$threshold"
