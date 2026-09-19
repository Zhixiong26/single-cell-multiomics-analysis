# Scanpy and Harmony

The portable implementation accepts 10x directories, ZIP archives containing exactly one matrix, and 10x H5. It performs per-sample QC, optional Scrublet, concatenation, normalization, HVG, PCA, optional Harmony, neighbours, UMAP, Leiden, markers, and guarded annotation.

RNA cell IDs use explicit `cell_id_prefix`. Preserve raw counts in `layers['counts']` and log-normalized full-gene values in `raw`. An annotation profile is valid only when its analysis signature and expected cluster set match exactly. Otherwise export an annotation template and `requires_review`; do not publish final cell-type claims.

`Notebooks/scanpy_workflow.ipynb` is the implementation, not a copy of one: the DAG's scanpy task runs `run_scanpy_notebook.py`, which writes a parameter sidecar and executes that same notebook through a kernel, so there is one code path rather than a notebook and a batch script that must be kept in step. `Scripts/Scanpy/README.md` documents the full ten-step notebook flow, the verification loop, and the completion standard; `Scripts/Scanpy/Report.md` is the bilingual evidence ledger. `Notebooks/example_ipf_scanpy.ipynb` is a shipped worked example that the DAG never executes — including its reviewed cluster-to-cell-type mapping, which is valid only for that example's exact cluster set and must not be transplanted.
