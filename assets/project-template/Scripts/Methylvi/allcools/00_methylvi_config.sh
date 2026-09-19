#!/usr/bin/env bash
# .cov -> ALLCools mCG 5-kb clustering -> MethylVI configuration.
# Export a variable before sourcing this file to override its default.

export SCMO_PROJECT_DIR="${SCMO_PROJECT_DIR:-${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}}"
export SCMO_METHSCAN_RUN_DIR="${SCMO_METHSCAN_RUN_DIR:?set SCMO_METHSCAN_RUN_DIR}"
export SCMO_INPUT_MANIFEST="${SCMO_INPUT_MANIFEST:-${SCMO_METHSCAN_RUN_DIR}/00_scanpy_selected/input_manifest.tsv}"
export SCMO_FILTERED_CELL_IDS="${SCMO_FILTERED_CELL_IDS:-${SCMO_METHSCAN_RUN_DIR}/03_filtered/column_header.txt}"
export SCMO_COV_DIR="${SCMO_COV_DIR:-${SCMO_PROJECT_DIR}/Data/COV}"
# Coverage and a pre-existing ALLC directory are fallback inputs only. The
# maintained route resolves original indexed ALLCs from the MethSCAn manifest.
export SCMO_EXISTING_ALLC_DIR="${SCMO_EXISTING_ALLC_DIR:?set SCMO_EXISTING_ALLC_DIR}"
export SCMO_ANNOTATION="${SCMO_ANNOTATION:?set SCMO_ANNOTATION}"
export SCMO_CHROM_SIZES="${SCMO_CHROM_SIZES:?set SCMO_CHROM_SIZES}"
export SCMO_BLACKLIST="${SCMO_BLACKLIST:?set SCMO_BLACKLIST}"
export SCMO_BLACKLIST_MD5="${SCMO_BLACKLIST_MD5:-}"
export SCMO_BLACKLIST_FRACTION="${SCMO_BLACKLIST_FRACTION:-0.2}"
# Feature selection changes the retained 5-kb features and therefore the
# MethylVI input, so each blacklist/fraction variant needs its own result root.
export SCMO_MVI_ROOT="${SCMO_MVI_ROOT:-${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi/allcools}"
export SCMO_ALLCOOLS_ROOT="${SCMO_ALLCOOLS_ROOT:-${SCMO_MVI_ROOT}/allcools_features}"
export SCMO_ALLC_DIR="${SCMO_ALLC_DIR:-${SCMO_ALLCOOLS_ROOT}/input_allc}"
export SCMO_ALLC_TABLE="${SCMO_ALLC_TABLE:-${SCMO_ALLCOOLS_ROOT}/selected_cells.allc.tsv}"
export SCMO_MCDS="${SCMO_MCDS:-${SCMO_ALLCOOLS_ROOT}/mcg.mcds}"
export SCMO_ALLCOOLS_H5AD="${SCMO_ALLCOOLS_H5AD:-${SCMO_ALLCOOLS_ROOT}/clustered.h5ad}"
export SCMO_MVI_INPUT="${SCMO_MVI_INPUT:-${SCMO_MVI_ROOT}/methylvi_input.h5mu}"
export SCMO_MVI_RESULTS="${SCMO_MVI_RESULTS:-${SCMO_MVI_ROOT}/results}"
# ALLCools and MethylVI use separate environments to avoid Python conflicts.
export SCMO_ALLCOOLS_ENV="${SCMO_ALLCOOLS_ENV:?set SCMO_ALLCOOLS_ENV}"
export SCMO_METHYLVI_ENV="${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}"
# Set this per route. Projects whose annotation table uses a different label
# column (for example a manually curated one) must set it explicitly.
export SCMO_CELLTYPE_KEY="${SCMO_CELLTYPE_KEY:-cell_type}"
# The integration batch key is project-defined; it must match the sample manifest.
export SCMO_BATCH_KEY="${SCMO_BATCH_KEY:-sample_id}"
export SCMO_BIN_SIZE="${SCMO_BIN_SIZE:-5000}"
export SCMO_MC_CONTEXT="${SCMO_MC_CONTEXT:-CGN}"
export SCMO_HYPO_SCORE_CUTOFF="${SCMO_HYPO_SCORE_CUTOFF:-0.9}"
export SCMO_BINARIZE_CUTOFF="${SCMO_BINARIZE_CUTOFF:-0.95}"
# The 30k branch is built first; the nested 10k input is derived from its
# deterministic prevalence ranking without rescanning ALLCs.
export SCMO_TARGET_FEATURES="${SCMO_TARGET_FEATURES:-30000}"
export SCMO_FEATURE_TARGETS="${SCMO_FEATURE_TARGETS:-10000 30000}"
export SCMO_LSI_COMPONENTS="${SCMO_LSI_COMPONENTS:-100}"
export SCMO_LSI_P_CUTOFF="${SCMO_LSI_P_CUTOFF:-0.1}"
export SCMO_ALLCOOLS_NEIGHBORS="${SCMO_ALLCOOLS_NEIGHBORS:-25}"
export SCMO_ALLCOOLS_LEIDEN_RESOLUTION="${SCMO_ALLCOOLS_LEIDEN_RESOLUTION:-1.0}"
export SCMO_CONSENSUS_LEIDEN_REPEATS="${SCMO_CONSENSUS_LEIDEN_REPEATS:-500}"
export SCMO_CONSENSUS_LEIDEN_RESOLUTION="${SCMO_CONSENSUS_LEIDEN_RESOLUTION:-0.5}"
# Derived from the filtered cell list at plan time, never a fixed constant.
# Cells without a manual label remain in the unsupervised analysis and receive
# an explicit Unknown label rather than being dropped.
export SCMO_EXPECTED_CELLS="${SCMO_EXPECTED_CELLS:?set dynamically derived SCMO_EXPECTED_CELLS}"
export SCMO_INCLUDE_UNANNOTATED="${SCMO_INCLUDE_UNANNOTATED:-1}"
export SCMO_MAX_CELLS="${SCMO_MAX_CELLS:-0}"
export SCMO_THREADS="${SCMO_THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
export SCMO_SEED="${SCMO_SEED:-0}"
export SCMO_EPOCHS="${SCMO_EPOCHS:-500}"
export SCMO_BATCH_SIZE="${SCMO_BATCH_SIZE:-32}"
# Supervised UMAP is calculated from the fixed MethylVI latent representation;
# Unknown labels remain unsupervised.
export SCMO_SUPERVISED_TARGET_KEY="${SCMO_SUPERVISED_TARGET_KEY:-${SCMO_CELLTYPE_KEY}}"
export SCMO_SUPERVISED_TARGET_WEIGHTS="${SCMO_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export SCMO_SUPERVISED_NEIGHBORS="${SCMO_SUPERVISED_NEIGHBORS:-15}"
export SCMO_SUPERVISED_MIN_DIST="${SCMO_SUPERVISED_MIN_DIST:-0.5}"
