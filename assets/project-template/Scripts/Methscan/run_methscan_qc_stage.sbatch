#!/usr/bin/env bash
# Run one serial MethSCAn QC stage with a small dedicated allocation.
# Usage: sbatch run_methscan_qc_stage.sbatch OUTPUT_DIR prepare|filter|smooth

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 2 )); then echo "Usage: sbatch $0 OUTPUT_DIR prepare|filter|smooth" >&2; exit 2; fi
output_dir=$1
stage=$2
case "$stage" in prepare|filter|smooth) ;; *) echo "Unknown QC stage: $stage" >&2; exit 3;; esac
[[ -f "$output_dir/common.COMPLETE" ]] || { echo "Common stage is not complete" >&2; exit 4; }
status_file="$output_dir/${stage}.status.tsv"
printf 'stage\tstatus\ttimestamp\n' > "$status_file"
record_stage() { printf '%s\t%s\t%s\n' "$stage" "$1" "$(date -Is)" >> "$status_file"; }
record_stage running
case "$stage" in
  prepare) "$methscan_python" "$script_dir/03_prepare_methscan.py" --methscan "$methscan_exe" --manifest "$output_dir/01_cov/input_manifest.tsv" --chunksize "$prepare_chunksize" --output-dir "$output_dir/02_prepared"; marker=prepare.COMPLETE ;;
  filter) "$methscan_exe" filter --min-sites "$min_sites" --min-meth "$min_meth" --max-meth "$max_meth" "$output_dir/02_prepared" "$output_dir/03_filtered"; marker=filter.COMPLETE ;;
  smooth) "$methscan_exe" smooth --bandwidth "$smooth_bandwidth" "$output_dir/03_filtered"; marker=smooth.COMPLETE ;;
esac
touch "$output_dir/$marker"
record_stage complete
printf 'QC stage completed: %s (%s)\n' "$output_dir" "$stage"
