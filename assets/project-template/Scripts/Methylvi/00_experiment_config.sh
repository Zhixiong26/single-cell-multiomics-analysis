#!/usr/bin/env bash
# Shared configuration for the MethylVI 10k/30k experiment.
# Export a variable before sourcing this file to override its default.
export MVI_PROJECT_DIR="${MVI_PROJECT_DIR:-${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}}"
export MVI_EXPERIMENT_ROOT="${MVI_EXPERIMENT_ROOT:-${SCMO_MVI_EXPERIMENT_ROOT:-${SCMO_RESULT_DIR:?set SCMO_RESULT_DIR}/methylvi}}"
export MVI_FEATURE_TARGETS="${MVI_FEATURE_TARGETS:-${SCMO_FEATURE_TARGETS:-10000 30000}}"
export MVI_VMR_THRESHOLDS="${MVI_VMR_THRESHOLDS:-${SCMO_VMR_THRESHOLDS:-0.01 0.02 0.05}}"
export MVI_MAX_FEATURES="${MVI_MAX_FEATURES:-${SCMO_MAX_FEATURES:?set SCMO_MAX_FEATURES}}"
export MVI_METHYLVI_ENV="${MVI_METHYLVI_ENV:-${SCMO_METHYLVI_ENV:?set SCMO_METHYLVI_ENV}}"
export MVI_ALLCOOLS_ENV="${MVI_ALLCOOLS_ENV:-${SCMO_ALLCOOLS_ENV:?set SCMO_ALLCOOLS_ENV}}"

# Re-export the resolved values under the neutral names so that every submit
# wrapper and stage script can read one namespace. SCMO_MVI_EXPERIMENT_ROOT is
# the parent of both routes; SCMO_MVI_ROOT (defined by the allcools route
# config) is only its allcools/ subdirectory.
export SCMO_MVI_EXPERIMENT_ROOT="$MVI_EXPERIMENT_ROOT"
export SCMO_FEATURE_TARGETS="$MVI_FEATURE_TARGETS"
export SCMO_VMR_THRESHOLDS="$MVI_VMR_THRESHOLDS"
export SCMO_MAX_FEATURES="$MVI_MAX_FEATURES"
export SCMO_METHYLVI_ENV="$MVI_METHYLVI_ENV"
export SCMO_ALLCOOLS_ENV="$MVI_ALLCOOLS_ENV"
