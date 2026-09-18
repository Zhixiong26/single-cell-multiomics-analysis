# ALLCools and MethylVI

Portable ALLCools-bin, MethSCAn-VMR, and VMR plus all-unique pooled-DMR routes. MethylVI uses integer `mc/cov`, validates `mc <= cov`, and exports ordinary/supervised UMAP plus methylation QC.

通用 ALLCools、VMR 和 VMR+全部 unique pooled-DMR 路线；训练使用整数 `mc/cov` 并输出完整 QC。

Production execution is available only through the project DAG (`tools/plan_workflow.py` and `tools/submit_workflow.py`); stage-local Slurm submission is intentionally unsupported.

正式任务只能通过项目 DAG（`tools/plan_workflow.py` 与 `tools/submit_workflow.py`）执行，不支持绕过证据链的阶段级 Slurm 直提。
