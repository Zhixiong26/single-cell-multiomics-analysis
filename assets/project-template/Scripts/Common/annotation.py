#!/usr/bin/env python3
"""Pure annotation-profile guard used by analysis stages and tests."""


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
