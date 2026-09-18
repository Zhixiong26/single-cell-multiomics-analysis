# ALLCools and MethylVI

ALLCools may rank blacklist-filtered bins for feature selection. MethylVI must use integer per-cell/per-feature `mc` and `cov`, not hypo-score.

Require unique cells/features, integral nonnegative counts, `mc <= cov`, explicit batch metadata, and a source manifest. Derive expected cells from the filtered-cell authority; zero is never a runtime default. Context, references, thresholds, feature targets, epochs, and batch size come from `analysis.yaml`. Nested feature targets come deterministically from one ranked maximum-feature input.

VMRs and ALLC manifests must share a completed MethSCAn signature. VMR+DMR retains selected VMRs and appends only merged pooled DMRs with zero genomic overlap.

Export latent embedding, ordinary UMAP/Leiden, supervised sensitivity views, methylation QC, summaries, and completion markers under the run-specific result root. Ordinary UMAP remains the primary latent-structure view.
