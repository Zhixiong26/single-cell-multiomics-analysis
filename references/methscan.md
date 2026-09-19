# MethSCAn, VMR, and DMR

The packaged workflow discovers indexed ALLCs independently for each `samples.tsv` row, applies the configured cell-ID rule and methylation context, then runs coverage conversion, prepare, filter, smooth, VMR scan/matrix, and Scanpy representation. Samples, thresholds, context, and chromosome order come from configuration and `chrom_sizes`; unmatched files are errors or explicit QC records, never silent skips.

Cell-type DMR requires `review_status: approved`, rejects placeholder labels, and enforces a minimum cell count. Per-sample and pooled-sample DMR are separate. Resume only comparisons whose identity, completion state, and format match.

`Scripts/Methscan/README.md` documents the stage sequence, output tree, and the `jobs × threads ≤ SLURM_CPUS_PER_TASK` constraint; `Scripts/Methscan/Report.md` is the bilingual evidence ledger. The `run_*.sbatch` wrappers and `submit_methscan_pipeline.sh` execute these stages standalone, sharing one dependency-chain implementation through `Scripts/Common/submit_helpers.sh`.

The heatmap branch may use raw-p fallback only for the known FDR division-by-zero failure. Top-N hypo-DMR heatmaps produce raw mean ratio, DMR-wise z-score, and display-compressed z-score figures. The MethylVI extension instead uses all pooled hypo-DMRs passing configured raw-p and absolute-difference thresholds, without Top-N truncation, and appends only zero-overlap DMRs to selected VMRs.
