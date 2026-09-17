# Scanpy and Harmony

The portable implementation accepts 10x directories, ZIP archives containing exactly one matrix, and 10x H5. It performs per-sample QC, optional Scrublet, concatenation, normalization, HVG, PCA, optional Harmony, neighbours, UMAP, Leiden, markers, and guarded annotation.

RNA cell IDs use explicit `cell_id_prefix`. Preserve raw counts in `layers['counts']` and log-normalized full-gene values in `raw`. An annotation profile is valid only when its analysis signature and expected cluster set match exactly. Otherwise export an annotation template and `requires_review`; do not publish final cell-type claims.

The clean notebook delegates to the same parameterized implementation as batch mode.
