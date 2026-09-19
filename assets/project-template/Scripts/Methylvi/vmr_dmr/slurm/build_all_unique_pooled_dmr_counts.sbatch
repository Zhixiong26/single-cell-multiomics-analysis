#!/usr/bin/env bash
# Count pooled-DMR coverage by reusing the VMR builder on the merged DMR set.
# Usage: sbatch [resource options] build_all_unique_pooled_dmr_counts.sbatch
# Resources and log paths are injected by the submitter, not declared here.
# The experiment configuration supplies the run root and the model ladder; it
# requires SCMO_MAX_FEATURES, and the reference paths below come from the caller.
# The route runner accepts a managed run ID, or the standalone acknowledgement
# this wrapper sets; see run.sh for the gate.
#SBATCH --job-name=scmo_pdmr_counts

set -Eeuo pipefail

# This wrapper is a standalone entry point: it acknowledges the non-DAG run
# that run.sh requires, so the stage can be explored, smoke tested, or
# recovered without a plan. Production submission still goes through the DAG.
export SCMO_STANDALONE_ACK=1
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi
source "$script_dir/vmr_dmr/00_vmr_dmr_config.sh"
bed="$SCMO_VMR_DMR_BED"
input="$SCMO_VMR_DMR_DMR_INPUT"
work="$SCMO_VMR_DMR_COUNT_ROWS"
# The cell set is identical across thresholds, so any branch supplies the
# selected-cell ALLC table; the first configured threshold is the canonical one.
allc="$SCMO_MVI_EXPERIMENT_ROOT/vmr/var_${SCMO_VMR_DMR_THRESHOLDS[0]}/shared/input/selected_cells.allc.tsv"
annotation="${SCMO_ANNOTATION:?set SCMO_ANNOTATION}"
[[ -f "$SCMO_VMR_DMR_SHARED/prepare.COMPLETE" ]] || { echo "DMR preparation incomplete" >&2; exit 3; }
features=$(wc -l < "$bed")
(( features >= 2 )) || { echo "Too few merged DMRs: $features" >&2; exit 4; }
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
# The VMR route runner sources its own configuration, so every path it requires
# is passed explicitly; the SCMO_* names are the ones the caller supplies.
VMR_FILTERED_BED="$bed" VMR_ALLC_TABLE="$allc" VMR_ANNOTATION="$annotation" \
VMR_COUNT_ROWS="$work" VMR_MVI_INPUT="$input" VMR_THREADS="${SLURM_CPUS_PER_TASK:-32}" \
VMR_MIN_COVERED_PERCENT=0 VMR_TARGET_FEATURES="$features" \
VMR_METHSCAN_RUN_DIR="${VMR_METHSCAN_RUN_DIR:-${SCMO_METHSCAN_RUN_DIR:?set SCMO_METHSCAN_RUN_DIR}}" \
VMR_EXISTING_ALLC_DIR="${VMR_EXISTING_ALLC_DIR:?set VMR_EXISTING_ALLC_DIR}" \
VMR_CHROM_SIZES="${VMR_CHROM_SIZES:-${SCMO_CHROM_SIZES:?set SCMO_CHROM_SIZES}}" \
VMR_BLACKLIST="${VMR_BLACKLIST:-${SCMO_BLACKLIST:?set SCMO_BLACKLIST}}" \
VMR_METHYLVI_ENV="${VMR_METHYLVI_ENV:-${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}}" \
  bash "$script_dir/vmr/run.sh" build
touch "$SCMO_VMR_DMR_SHARED/counts.COMPLETE"
