#!/usr/bin/env bash
# Portable defaults. Project-specific paths are injected by the generated task command.

project_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
methscan_exe=${SCMO_METHSCAN_EXE:?set SCMO_METHSCAN_EXE}
methscan_python=${SCMO_METHSCAN_PYTHON:?set SCMO_METHSCAN_PYTHON}
scanpy_python=${SCMO_SCANPY_PYTHON:?set SCMO_SCANPY_PYTHON}
rna_annotation_table=${SCMO_ANNOTATION:-}
rna_exclude_cell_type=${SCMO_EXCLUDE_CELL_TYPE:-NA}
allc_source=${SCMO_ALLC_SOURCE:-}
read -r -a samples <<< "${SCMO_SAMPLE_IDS:-}"

min_sites=${SCMO_MIN_SITES:-300000}
min_meth=${SCMO_MIN_METH_PERCENT:-50}
max_meth=${SCMO_MAX_METH_PERCENT:-100}
smooth_bandwidth=${SCMO_SMOOTH_BANDWIDTH:-1000}
scan_bandwidth=${SCMO_SCAN_BANDWIDTH:-2000}
scan_stepsize=${SCMO_SCAN_STEPSIZE:-100}
read -r -a scan_var_thresholds <<< "${SCMO_VMR_THRESHOLDS:-0.01 0.02 0.05}"
scan_min_cells=${SCMO_MIN_CELLS:-6}
methdiff_min_cells=${SCMO_MIN_CELLS:-6}
methdiff_bandwidth=${SCMO_DIFF_BANDWIDTH:-2000}
methdiff_stepsize=${SCMO_DIFF_STEPSIZE:-1000}
methdiff_threshold=${SCMO_DIFF_THRESHOLD:-0.02}
pooled_sample_label=${SCMO_POOLED_LABEL:-pooled_samples}
methdiff_jobs=${SCMO_METHDIFF_JOBS:-2}
methdiff_threads=${SCMO_METHDIFF_THREADS:-4}
hypo_raw_p=${SCMO_HYPO_RAW_P:-0.01}
hypo_min_abs_diff=${SCMO_HYPO_MIN_ABS_DIFF:-0.25}
hypo_top_per_cell_type=${SCMO_HYPO_TOP_PER_CELL_TYPE:-200}
hypo_matrix_workers=${SCMO_HYPO_MATRIX_WORKERS:-8}
hypo_zscore_min_observed_cells=${SCMO_ZSCORE_MIN_OBSERVED_CELLS:-30}
hypo_zscore_clip=${SCMO_ZSCORE_CLIP:-3}
hypo_compressed_colorbar_limit=${SCMO_COMPRESSED_COLORBAR_LIMIT:-1}
hypo_heatmap_dpi=${SCMO_HEATMAP_DPI:-300}
threads=${SLURM_CPUS_PER_TASK:-${SCMO_THREADS:-16}}
prepare_chunksize=${SCMO_PREPARE_CHUNKSIZE:-10000000}
min_free_gb=${SCMO_MIN_FREE_GB:-10}
cov_conversion_workers=${SCMO_COV_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}
cov_compresslevel=${SCMO_COV_COMPRESSLEVEL:-1}
vmr_min_cell_fraction=${SCMO_VMR_MIN_CELL_FRACTION:-0.05}
cell_min_regions=${SCMO_CELL_MIN_REGIONS:-100}
pca_components=${SCMO_PCA_COMPONENTS:-30}
neighbors=${SCMO_NEIGHBORS:-15}
leiden_resolution=${SCMO_LEIDEN_RESOLUTION:-0.8}
random_seed=${SCMO_SEED:-0}
