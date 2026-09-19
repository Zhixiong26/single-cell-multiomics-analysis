# ALLCools and MethylVI

ALLCools may rank blacklist-filtered bins for feature selection. MethylVI must use integer per-cell/per-feature `mc` and `cov`, not hypo-score.

Require unique cells/features, integral nonnegative counts, `mc <= cov`, explicit batch metadata, and a source manifest. Derive expected cells from the filtered-cell authority; zero is never a runtime default. Context, references, thresholds, feature targets, epochs, batch size, and the validation fraction come from `analysis.yaml`. Nested feature targets come deterministically from one ranked maximum-feature input. `methylvi.validation_fraction` reaches scvi as `validation_size`, and `train_size` is set to `1 - validation_fraction` alongside it: scvi defaults `train_size` to 0.9 on its own and rejects `train_size + validation_size > 1`, so passing only the validation size would fail for any fraction above 0.1.

VMRs and ALLC manifests must share a completed MethSCAn signature. VMR+DMR retains selected VMRs and appends only merged pooled DMRs with zero genomic overlap.

Export latent embedding, ordinary UMAP/Leiden, supervised sensitivity views, methylation QC, summaries, and completion markers under the run-specific result root. Ordinary UMAP remains the primary latent-structure view.

ALLCools reuses the common MethSCAn selection/prepare/filter chain because its expected-cell authority and COV inputs come from the filtered manifest. It does not require smooth or VMR scan unless a VMR route is also selected.

Three routes share one training core: `allcools` (5-kb bins), `vmr` (variance-thresholded VMR regions), and `vmr_dmr` (selected VMRs plus every unique pooled DMR). Each route owns a `00_*_config.sh` declaring its directory layout and defaults; `vmr_dmr` reuses the `vmr` feature inputs rather than rescanning ALLCs, so it depends on that route's completed signature. Route configs resolve the shared experiment values from the entry `00_experiment_config.sh` and re-export them under the neutral `SCMO_*` names, so a stage reads one namespace regardless of which route invoked it. `SCMO_MVI_EXPERIMENT_ROOT` is the parent of all routes; `SCMO_MVI_ROOT` is only its `allcools/` subdirectory.

The packaged `run.sh` runners and `slurm/*.sbatch` wrappers execute one stage without a plan, for exploration, smoke testing, and recovery. They resolve resources from scheduler profiles, verify the same completion markers the DAG does, and never substitute for a run's completion evidence. Production MethylVI work goes through the run DAG.
