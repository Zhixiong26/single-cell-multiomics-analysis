#!/usr/bin/env bash
# Run selection and ALLC-to-COV only. Later QC stages are separate jobs.
# Usage: sbatch [resource options] run_methscan_common.sbatch OUTPUT_DIR

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
[[ "$(readlink -m "$output_dir")" == "$(readlink -m "${SCMO_METHSCAN_RUN_DIR:?set SCMO_METHSCAN_RUN_DIR}")" ]] || { echo "Unexpected MethSCAn output root" >&2; exit 3; }
[[ ! -e "$output_dir" ]] || { echo "Output already exists: $output_dir" >&2; exit 4; }
[[ -x "$methscan_exe" && -x "$methscan_python" && -x "$scanpy_python" ]] || { echo "Configured executable missing" >&2; exit 5; }
[[ -r "${SCMO_SAMPLES_TSV:?set SCMO_SAMPLES_TSV}" ]] || { echo "Input path unavailable" >&2; exit 6; }
mkdir -p "$(dirname "$output_dir")"
available_kb=$(df -Pk "$(dirname "$output_dir")" | awk 'NR==2 {print $4}')
required_kb=$((min_free_gb * 1024 * 1024))
(( available_kb >= required_kb )) || { echo "Insufficient free space" >&2; exit 7; }
mkdir -p "$output_dir"
status_file="$output_dir/common_status.tsv"
printf 'stage\tstatus\ttimestamp\n' > "$status_file"
{
  printf 'key\tvalue\n'
  printf 'samples_tsv\t%s\n' "$SCMO_SAMPLES_TSV"
  printf 'samples_tsv_sha256\t%s\n' "$(sha256sum "$SCMO_SAMPLES_TSV" | awk '{print $1}')"
  printf 'sample_ids\t%s\n' "$sample_ids"
  printf 'allc_sources_json\t%s\n' "$allc_sources"
  printf 'rna_annotation_table\t%s\n' "$rna_annotation_table"
  printf 'scan_var_thresholds\t%s\n' "${scan_var_thresholds[*]}"
  printf 'methdiff_scope\tper_sample\n'
  printf 'methdiff_min_cells\t%s\n' "$methdiff_min_cells"
  printf 'hypo_raw_p_strictly_less_than\t%s\n' "$hypo_raw_p"
  printf 'hypo_min_abs_diff\t%s\n' "$hypo_min_abs_diff"
  printf 'hypo_top_per_cell_type\t%s\n' "$hypo_top_per_cell_type"
  printf 'hypo_zscore_clip\t%s\n' "$hypo_zscore_clip"
  printf 'hypo_compressed_colorbar_limit\t%s\n' "$hypo_compressed_colorbar_limit"
} > "$output_dir/effective_configuration.tsv"
record_stage() { printf '%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" >> "$status_file"; }
current_stage=initialization
on_error() { local code=$?; trap - ERR; record_stage "$current_stage" "failed(exit=$code)"; exit "$code"; }
trap on_error ERR
run_stage() { current_stage=$1; shift; record_stage "$current_stage" running; "$methscan_python" "$script_dir/06_run_with_resources.py" --output "$output_dir/${current_stage}.resources.json" -- "$@"; record_stage "$current_stage" complete; }
{
  "$methscan_exe" --version; "$methscan_python" --version
  "$scanpy_python" -c 'import scanpy; print("scanpy=" + scanpy.__version__)'
} > "$output_dir/software_versions.txt"
selection_args=(--samples-tsv "$SCMO_SAMPLES_TSV" --exclude-cell-type "$rna_exclude_cell_type" --output-dir "$output_dir/00_scanpy_selected")
if [[ -n "$rna_annotation_table" ]]; then selection_args+=(--scanpy-annotation "$rna_annotation_table"); fi
if [[ "${SCMO_ANNOTATION_APPROVED:-0}" == 1 ]]; then selection_args+=(--annotation-approved); fi
run_stage 00_scanpy_selection "$methscan_python" "$script_dir/01_select_scanpy_cells.py" "${selection_args[@]}"
run_stage 01_convert_cov "$methscan_python" "$script_dir/02_convert_allc_to_cov.py" --manifest "$output_dir/00_scanpy_selected/input_manifest.tsv" --output-dir "$output_dir/01_cov" --workers "${SLURM_CPUS_PER_TASK:-$cov_conversion_workers}" --compresslevel "$cov_compresslevel" --mc-context "$SCMO_MC_CONTEXT"
touch "$output_dir/common.COMPLETE"
record_stage common workflow_complete
printf 'Common MethSCAn stages completed: %s\n' "$output_dir"
