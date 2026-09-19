#!/usr/bin/env bash
# Submit the pooled all-sample pairwise cell-type DMR analysis only.
# Standalone opt-in branch for a run whose common and smooth stages completed.
# Usage: submit_methscan_pooled_dmr.sh OUTPUT_DIR
set -Eeuo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SCMO_PROJECT_ROOT=${SCMO_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
source "$script_dir/00_methscan_config.sh"
source "$script_dir/../Common/submit_helpers.sh"
if (( $# != 1 )); then echo "Usage: $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
diff_job=$(scmo_submit dmr "" "$script_dir/run_methscan_pooled_methdiff.sbatch" "$output_dir")
printf 'pooled_methdiff=%s\n' "$diff_job"
