#!/usr/bin/env bash
# VMR + all-unique-pooled-DMR route configuration.
#
# This route reuses the VMR feature inputs built by the vmr route and augments
# them with the union of every unique hypo-DMR from the pooled MethSCAn run, so
# it depends on that run's pairwise summary rather than on ALLC rescanning.
#
# This file declares the route's directory layout and its defaults. It requires
# only what the route cannot run without; each stage requires the specific inputs
# it consumes, so a summary stage is not blocked by an unrelated missing path.
# Export a variable before sourcing this file to override its default.

export SCMO_VMR_DMR_SCRIPT_DIR="${SCMO_VMR_DMR_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$SCMO_VMR_DMR_SCRIPT_DIR/../00_experiment_config.sh"

# Route root. The shared/ subtree holds the merged DMR inputs and their counts;
# each var_/features_ leaf below it holds one combined training input.
export SCMO_VMR_DMR_ROOT="${SCMO_VMR_DMR_ROOT:-${SCMO_MVI_EXPERIMENT_ROOT}/vmr_dmr}"
export SCMO_VMR_DMR_SHARED="${SCMO_VMR_DMR_SHARED:-${SCMO_VMR_DMR_ROOT}/shared/dmr}"
export SCMO_VMR_DMR_BED="${SCMO_VMR_DMR_BED:-${SCMO_VMR_DMR_SHARED}/all_unique_hypo_DMRs.merged.bed}"
export SCMO_VMR_DMR_DMR_INPUT="${SCMO_VMR_DMR_DMR_INPUT:-${SCMO_VMR_DMR_SHARED}/all_unique_hypo_DMRs.input.h5mu}"
export SCMO_VMR_DMR_COUNT_ROWS="${SCMO_VMR_DMR_COUNT_ROWS:-${SCMO_VMR_DMR_SHARED}/count_rows}"

# leaf_root THRESHOLD FEATURES -> directory holding one input.h5mu and results/.
leaf_root() {
  printf '%s/var_%s/features_%s\n' "$SCMO_VMR_DMR_ROOT" "$1" "$2"
}
# vmr_input THRESHOLD FEATURES -> upstream VMR feature input this route extends.
# That input lives under the vmr route root, a sibling of this route's root.
vmr_input() {
  printf '%s/vmr/var_%s/features_%s/input.h5mu\n' "$SCMO_MVI_EXPERIMENT_ROOT" "$1" "$2"
}

# The model ladder is validated against these lists instead of a hardcoded case
# statement, so a new threshold or feature target needs no script edit.
read -r -a SCMO_VMR_DMR_THRESHOLDS <<< "$SCMO_VMR_THRESHOLDS"
read -r -a SCMO_VMR_DMR_FEATURE_TARGETS <<< "$SCMO_FEATURE_TARGETS"

# Hypo-DMR selection thresholds. SCMO_HYPO_RAW_P and SCMO_HYPO_MIN_ABS_DIFF are
# the established names: the run DAG injects them from analysis.methscan.dmr, and
# the MethSCAn hypo-heatmap stages read the same values. They must stay identical
# between the prepare and count stages, or the merged BED will not match the
# counts built from it.
export SCMO_HYPO_RAW_P="${SCMO_HYPO_RAW_P:-0.01}"
export SCMO_HYPO_MIN_ABS_DIFF="${SCMO_HYPO_MIN_ABS_DIFF:-0.25}"
export SCMO_BLACKLIST_FRACTION="${SCMO_BLACKLIST_FRACTION:-0.2}"

# Exported for the shared MethylVI training and UMAP code.
export SCMO_BATCH_KEY="${SCMO_BATCH_KEY:-sample_id}"
export SCMO_CELLTYPE_KEY="${SCMO_CELLTYPE_KEY:-cell_type}"
export SCMO_SUPERVISED_TARGET_KEY="${SCMO_SUPERVISED_TARGET_KEY:-$SCMO_CELLTYPE_KEY}"
export SCMO_SUPERVISED_TARGET_WEIGHTS="${SCMO_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export SCMO_SUPERVISED_NEIGHBORS="${SCMO_SUPERVISED_NEIGHBORS:-15}"
export SCMO_SUPERVISED_MIN_DIST="${SCMO_SUPERVISED_MIN_DIST:-0.5}"
