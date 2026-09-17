#!/usr/bin/env python3
"""Validate all three threshold branches and write machine-readable summaries."""
import argparse, csv, gzip, json
from pathlib import Path

MATRIX_FILES = ("methylated_sites.csv.gz", "total_sites.csv.gz", "methylation_fractions.csv.gz", "mean_shrunken_residuals.csv.gz")

def tsv(path):
    with path.open(newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def lines(path):
    with path.open() as h: return [x.strip() for x in h if x.strip()]
def ids(path, delimiter=",", gz=False):
    op = gzip.open if gz else open
    with op(str(path), "rt", newline="") as h:
        r = csv.reader(h, delimiter=delimiter); next(r, None); return [x[0] for x in r if x]
def required(path):
    if not path.is_file() or path.stat().st_size == 0: raise FileNotFoundError(path)
def unique(name, values):
    if len(values) != len(set(values)): raise ValueError("Duplicate IDs in " + name)
def vmr_count(path):
    with path.open() as h: return sum(bool(x.strip()) and not x.startswith("#") for x in h)

def main():
    p = argparse.ArgumentParser(); p.add_argument("--run-dir", type=Path, required=True); p.add_argument("--threshold", action="append", required=True); a=p.parse_args(); root=a.run_dir.resolve()
    common = {
        "selection": root/"00_scanpy_selected/input_manifest.tsv",
        "selection_summary": root/"00_scanpy_selected/scanpy_selection_summary.json",
        "cov_manifest": root/"01_cov/input_manifest.tsv", "cpg_qc": root/"01_cov/cell_cpg_qc.tsv",
        "conversion": root/"01_cov/conversion_summary.json", "prepared_stats": root/"02_prepared/cell_stats.csv",
        "filtered_stats": root/"03_filtered/cell_stats.csv", "prepared_header": root/"02_prepared/column_header.txt",
        "filtered_header": root/"03_filtered/column_header.txt", "id_check": root/"02_prepared/cell_id_check.json",
        "methdiff_summary": root/"07_methdiff/pairwise_summary.json",
        "hypo_heatmap_summary": root/"hypo_dmr_heatmap_summary.json",
    }
    for path in common.values(): required(path)
    methdiff = json.loads(common["methdiff_summary"].read_text())
    if methdiff.get("status") not in {"complete", "complete_with_fallback_needed"}:
        raise ValueError("Pairwise meth-diff summary has unresolved hard failures")
    hypo_heatmaps = json.loads(common["hypo_heatmap_summary"].read_text())
    if hypo_heatmaps.get("status") != "complete":
        raise ValueError("Hypo-DMR heatmap summary is not complete")
    for sample in hypo_heatmaps.get("plots", {}).get("samples", []):
        if sample.get("status") == "complete":
            for key in ("mean_ratio_plot", "zscore_plot", "zscore_compressed_colorbar_plot"):
                required(Path(sample[key]))
    selected=tsv(common["selection"]); selected_ids=[r["cell_id"] for r in selected]; unique("selection", selected_ids)
    if any(not (r.get("rna_cell_type") or "").strip() or (r.get("rna_cell_type") or "").strip()=="NA" for r in selected): raise ValueError("Selected manifest contains empty/NA cell types")
    if json.loads(common["id_check"].read_text()).get("status") != "pass": raise ValueError("prepare cell-ID check failed")
    prepared=lines(common["prepared_header"]); prepared_stats=[r["cell_name"] for r in csv.DictReader(common["prepared_stats"].open(newline=""))]
    filtered=lines(common["filtered_header"]); filtered_stats=[r["cell_name"] for r in csv.DictReader(common["filtered_stats"].open(newline=""))]
    if prepared != sorted(selected_ids) or prepared_stats != prepared or filtered_stats != filtered or not set(filtered).issubset(prepared): raise ValueError("prepare/filter cell IDs are inconsistent")
    branches={}
    for threshold in a.threshold:
        label="var_"+threshold; scan=root/"04_scan"/label/"VMRs.bed"; matrix=root/"05_matrix"/label; sp=root/"06_scanpy"/label
        required(scan)
        for name in MATRIX_FILES: required(matrix/name)
        for name in ("tables/cell_embedding.tsv","tables/cell_missingness_qc.tsv","tables/vmr_missingness_qc.tsv","objects/methscan_vmr.h5ad","run_parameters.json","figures/umap_leiden.png","figures/umap_rna_cell_type.png"): required(sp/name)
        for name in MATRIX_FILES:
            mid=ids(matrix/name, gz=True); unique(label+"/"+name, mid)
            if mid != filtered: raise ValueError(label+" matrix IDs differ from filtered IDs")
        emb=tsv(sp/"tables/cell_embedding.tsv"); emb_ids=[r[next(iter(emb[0]))] for r in emb] if emb else []; unique(label+" embedding", emb_ids)
        params=json.loads((sp/"run_parameters.json").read_text())
        if int(params["input_cells"]) != len(filtered) or int(params["retained_cells"]) != len(emb_ids) or not set(emb_ids).issubset(filtered): raise ValueError(label+" Scanpy counts/IDs inconsistent")
        branches[label]={"threshold":float(threshold),"vmrs":vmr_count(scan),"matrix_cells":len(filtered),"scanpy_input_cells":int(params["input_cells"]),"scanpy_cells":int(params["retained_cells"]),"scanpy_vmrs":int(params["retained_regions"]),"vmr_bed":str(scan),"matrix_dir":str(matrix),"scanpy_dir":str(sp)}
    selection=json.loads(common["selection_summary"].read_text()); filtered_n=len(filtered); selected_n=len(selected); input_n=int(selection["discovered_allc_cells"])
    summary={"status":"complete","run_dir":str(root),"canonical_cell_id":"<sample_id>_<17bp_barcode>","input_cells":input_n,"scanpy_selected_cells":selected_n,"prepared_cells":len(prepared),"filtered_cells":filtered_n,"thresholds":[float(x) for x in a.threshold],"branches":branches,"methdiff":methdiff,"hypo_dmr_heatmaps":hypo_heatmaps,"filter_retention_fraction":filtered_n/selected_n,"overall_allc_retention_fraction":filtered_n/input_n}
    (root/"run_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    with (root/"run_summary.tsv").open("w") as h:
        h.write("threshold\tvmrs\tmatrix_cells\tscanpy_input_cells\tscanpy_cells\tscanpy_vmrs\n")
        for b in branches.values(): h.write("{threshold}\t{vmrs}\t{matrix_cells}\t{scanpy_input_cells}\t{scanpy_cells}\t{scanpy_vmrs}\n".format(**b))
    print(json.dumps(summary,indent=2,sort_keys=True))
if __name__ == "__main__": main()
