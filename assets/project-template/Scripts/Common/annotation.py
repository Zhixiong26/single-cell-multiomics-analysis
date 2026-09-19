#!/usr/bin/env python3
"""Pure annotation-profile guard used by analysis stages and tests."""


# Labels that name the absence of a cell type rather than a cell type. The analysis stages each
# carry their own inline copy of the narrower subset they skip when reading metadata
# (Scripts/Methscan/01_select_scanpy_cells.py, 07_methdiff_celltype.py, and the supervised UMAP in
# Scripts/Methylvi/shared/05_plot_supervised_umap.py); this is the superset a human review must not
# be able to record as a cell type, because pooling or dropping those cells is not what the label
# says. Compared case-insensitively.
PLACEHOLDER_LABELS = frozenset({
    "", "na", "nan", "none", "n/a", "unknown", "unannotated", "unassigned",
    "requires_review", "requires review",
})


def evaluate_annotation_profile(observed_clusters, analysis_signature, profile, unassigned="Unassigned"):
    observed = sorted(str(item) for item in observed_clusters)
    template = {
        "analysis_signature": analysis_signature,
        "expected_clusters": observed,
        "annotations": {
            cluster: {"cell_type": unassigned, "confidence": "unreviewed", "evidence": ""}
            for cluster in observed
        },
    }
    status = {"status": "requires_review", "analysis_signature": analysis_signature, "observed_clusters": observed}
    if profile is None:
        status["reason"] = "annotation_profile_absent"
        return {}, status, template
    expected = sorted(str(item) for item in profile.get("expected_clusters", []))
    if profile.get("analysis_signature") != analysis_signature:
        status["reason"] = "analysis_signature_mismatch"
        return {}, status, template
    if expected != observed:
        status["reason"] = "cluster_set_mismatch"
        return {}, status, template
    annotations = {str(key): value for key, value in (profile.get("annotations") or {}).items()}
    if set(annotations) != set(observed) or any(not item.get("cell_type") for item in annotations.values()):
        status["reason"] = "annotation_mapping_incomplete"
        return {}, status, template
    status["status"] = "reviewed"
    return annotations, status, template
