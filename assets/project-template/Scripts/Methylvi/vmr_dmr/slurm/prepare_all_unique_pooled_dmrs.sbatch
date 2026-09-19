#!/usr/bin/env bash
# Select all unique pooled hypo-DMRs and merge overlapping DMR intervals.
# Usage: sbatch [resource options] prepare_all_unique_pooled_dmrs.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The experiment configuration supplies the run root and the model ladder; it
# requires SCMO_MAX_FEATURES, and the reference paths below come from the caller.
#SBATCH --job-name=scmo_pdmr_prepare

set -Eeuo pipefail
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi/vmr_dmr
source "$script_dir/00_vmr_dmr_config.sh"
output="$SCMO_VMR_DMR_SHARED"
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
"${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}/bin/python" "$script_dir/01_prepare_all_unique_pooled_dmrs.py" \
  --pairwise-summary "${SCMO_METHSCAN_RUN_DIR:?set SCMO_METHSCAN_RUN_DIR}/pooled/07_methdiff/pairwise_summary.tsv" \
  --blacklist "${SCMO_BLACKLIST:?set SCMO_BLACKLIST}" \
  --chrom-sizes "${SCMO_CHROM_SIZES:?set SCMO_CHROM_SIZES}" \
  --blacklist-fraction "$SCMO_BLACKLIST_FRACTION" \
  --raw-p "$SCMO_HYPO_RAW_P" --min-abs-diff "$SCMO_HYPO_MIN_ABS_DIFF" \
  --sort-threads "${SLURM_CPUS_PER_TASK:-8}" --output-dir "$output"
touch "$output/prepare.COMPLETE"
