# MethSCAn, VMR, and DMR

The packaged workflow selects indexed ALLCs from explicit metadata, converts selected CG records to Bismark coverage, then runs prepare, filter, smooth, VMR scan/matrix, and Scanpy representation. Samples and thresholds come from configuration.

Cell-type DMR requires reviewed labels and a minimum cell count. Per-sample and pooled-sample DMR are separate. Resume only comparisons whose identity, completion state, and format match.

The heatmap branch may use raw-p fallback only for the known FDR division-by-zero failure. Top-N hypo-DMR heatmaps produce raw mean ratio, DMR-wise z-score, and display-compressed z-score figures. The MethylVI extension instead uses all pooled hypo-DMRs passing configured raw-p and absolute-difference thresholds, without Top-N truncation, and appends only zero-overlap DMRs to selected VMRs.
