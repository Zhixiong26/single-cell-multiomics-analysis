# Documentation and evidence

Generated projects and major stages have bilingual README and Report files. README holds stable purpose, inputs, outputs, environments, parameters, order, entry points, and completion criteria. Report holds observed commands, jobs, resources, signatures, QC, failures, recovery, limitations, and unresolved decisions.

Keep dynamic state under `.workflow/runs/<run-id>/`. A stage is complete only when summaries, artifacts, and completion markers agree; Slurm `COMPLETED` alone is insufficient. Figures require adjacent tabular or JSON QC/manifests.

Do not rewrite README with runtime state. `update_report.py` updates Report; `inspect_run.py` writes machine-readable status. Analysis artifacts remain under the corresponding `Results/runs/<run-id>/` root.

Environment discovery, creation, verification, and failures are recorded under `.workflow/environment-bootstrap/` and summarized in `Scripts/Environment/Report.md`. Never report a profile as ready from package installation alone; its profile verification command must pass.
