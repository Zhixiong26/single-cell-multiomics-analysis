#!/usr/bin/env bash
# Independent VMR-region MethylVI configuration.
# Export a variable before sourcing this file to override its default.

export VMR_PROJECT_DIR="${VMR_PROJECT_DIR:-${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}}"
export VMR_SCRIPT_DIR="${VMR_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# The MethSCAn run supplies selected original ALLC paths; the filtered cell
# header below restricts them to the post-filter cell set.
export VMR_METHSCAN_RUN_DIR="${VMR_METHSCAN_RUN_DIR:?set VMR_METHSCAN_RUN_DIR}"
export VMR_METHSCAN_VARIANCE="${VMR_METHSCAN_VARIANCE:-0.01}"
export VMR_INPUT_MANIFEST="${VMR_INPUT_MANIFEST:-${VMR_METHSCAN_RUN_DIR}/00_scanpy_selected/input_manifest.tsv}"
export VMR_FILTERED_CELL_IDS="${VMR_FILTERED_CELL_IDS:-${VMR_METHSCAN_RUN_DIR}/03_filtered/column_header.txt}"
# Set this to the selected MethSCAn branch, for example:
# Results/Methscan/<run>/04_scan/var_0.01/VMRs.bed
export VMR_SOURCE_BED="${VMR_SOURCE_BED:-${VMR_METHSCAN_RUN_DIR}/04_scan/var_${VMR_METHSCAN_VARIANCE}/VMRs.bed}"
export VMR_COV_DIR="${VMR_COV_DIR:-${VMR_PROJECT_DIR}/Data/COV}"
export VMR_EXISTING_ALLC_DIR="${VMR_EXISTING_ALLC_DIR:?set VMR_EXISTING_ALLC_DIR}"
export VMR_CHROM_SIZES="${VMR_CHROM_SIZES:?set VMR_CHROM_SIZES}"
export VMR_BLACKLIST="${VMR_BLACKLIST:?set VMR_BLACKLIST}"
export VMR_BLACKLIST_MD5="${VMR_BLACKLIST_MD5:-}"
export VMR_BLACKLIST_FRACTION="${VMR_BLACKLIST_FRACTION:-0.2}"
export VMR_ANNOTATION="${VMR_ANNOTATION:?set VMR_ANNOTATION}"
export VMR_RESULTS_ROOT="${VMR_RESULTS_ROOT:-${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi/vmr/var_${VMR_METHSCAN_VARIANCE}}"
export VMR_INPUT_DIR="${VMR_INPUT_DIR:-${VMR_RESULTS_ROOT}/input}"
export VMR_FILTERED_BED="${VMR_FILTERED_BED:-${VMR_INPUT_DIR}/vmr.filtered.bed}"
export VMR_ALLC_TABLE="${VMR_ALLC_TABLE:-${VMR_INPUT_DIR}/selected_cells.allc.tsv}"
export VMR_MVI_INPUT="${VMR_MVI_INPUT:-${VMR_RESULTS_ROOT}/methylvi_input.h5mu}"
export VMR_COUNT_ROWS="${VMR_COUNT_ROWS:-${VMR_RESULTS_ROOT}/count_rows}"
export VMR_MVI_RESULTS="${VMR_MVI_RESULTS:-${VMR_RESULTS_ROOT}/results}"
export VMR_METHYLVI_ENV="${VMR_METHYLVI_ENV:-${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}}"
# Set a positive value only for the legacy coverage-directory fallback.
export VMR_EXPECTED_CELLS="${VMR_EXPECTED_CELLS:-0}"
# Optional coverage floor. The effective cutoff is calculated from the current
# cells as the highest covered-cell threshold that still permits the 30k target.
export VMR_MIN_COVERED_PERCENT="${VMR_MIN_COVERED_PERCENT:-0}"
export VMR_TARGET_FEATURES="${VMR_TARGET_FEATURES:-30000}"
export VMR_FEATURE_TARGETS="${VMR_FEATURE_TARGETS:-10000 30000}"
export VMR_THREADS="${VMR_THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
export VMR_SEED="${VMR_SEED:-0}"
export VMR_EPOCHS="${VMR_EPOCHS:-500}"
export VMR_BATCH_SIZE="${VMR_BATCH_SIZE:-32}"
export VMR_BATCH_KEY="${VMR_BATCH_KEY:-sample_id}"
export VMR_CELLTYPE_KEY="${VMR_CELLTYPE_KEY:-cell_type}"
# Compatibility variables consumed by the shared MethylVI training/UMAP code.
# Guarded with `:-` so that an explicitly exported SCMO_* value is not silently
# clobbered when this file is sourced after another route's configuration.
export SCMO_MVI_INPUT="${SCMO_MVI_INPUT:-$VMR_MVI_INPUT}"
export SCMO_MVI_RESULTS="${SCMO_MVI_RESULTS:-$VMR_MVI_RESULTS}"
export SCMO_THREADS="${SCMO_THREADS:-$VMR_THREADS}"
export SCMO_SEED="${SCMO_SEED:-$VMR_SEED}"
export SCMO_EPOCHS="${SCMO_EPOCHS:-$VMR_EPOCHS}"
export SCMO_BATCH_SIZE="${SCMO_BATCH_SIZE:-$VMR_BATCH_SIZE}"
export SCMO_BATCH_KEY="${SCMO_BATCH_KEY:-$VMR_BATCH_KEY}"
export SCMO_CELLTYPE_KEY="${SCMO_CELLTYPE_KEY:-$VMR_CELLTYPE_KEY}"
export SCMO_SUPERVISED_TARGET_KEY="${SCMO_SUPERVISED_TARGET_KEY:-$VMR_CELLTYPE_KEY}"
export SCMO_SUPERVISED_TARGET_WEIGHTS="${SCMO_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export SCMO_SUPERVISED_NEIGHBORS="${SCMO_SUPERVISED_NEIGHBORS:-15}"
export SCMO_SUPERVISED_MIN_DIST="${SCMO_SUPERVISED_MIN_DIST:-0.5}"
