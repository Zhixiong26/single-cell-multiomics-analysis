# Documentation and evidence

Generated projects and major stages have bilingual README and Report files. README holds stable purpose, inputs, outputs, environments, parameters, order, entry points, and completion criteria. Report holds observed commands, jobs, resources, signatures, QC, failures, recovery, limitations, and unresolved decisions.

Two levels exist, and they are owned differently. The root `README.md`/`Report.md` are the project's contract and run log: the README is written once by `init_project.py` from the intake and the resolved configuration, and the Report accumulates one concise record per run. The module documents under `Scripts/<Module>/` (and `Scripts/Methylvi/<route>/`) are static reference material ported from the template; no run updates them, and they must not be used to record project state.

The root README carries one section the module documents do not: the `Data/` entries and where each came from. It is generated with the project and frozen with the rest of the file, so it describes the inputs as configured rather than as they were later moved. The root Report repeats that table above the run-log region together with two things nothing else states — which bundled reference files generation filled in, and any sample whose input was configured as a path outside the project. Both are above the region and are therefore human text the tools never rewrite.

Keep dynamic state under `.workflow/runs/<run-id>/`. A stage is complete only when summaries, artifacts, and completion markers agree; Slurm `COMPLETED` alone is insufficient. Figures require adjacent tabular or JSON QC/manifests.

Never rewrite the root README with runtime state: it is frozen after generation, because a document that reports progress stops being trustworthy as a description of what the project is. Runtime facts belong in the run log inside `Report.md`. `update_report.py` rebuilds that log and `inspect_run.py` refreshes the inspected run's record as a side effect of computing the run's verdict.

Documentation writes never change a run's verdict. The exit code of `inspect_run.py` encodes the state of the run, not of a derived markdown file, so an unwritable Report is reported on stderr with the recovery command and leaves the exit code alone; `update_report.py` is the loud path and does exit non-zero. Analysis artifacts remain under the corresponding `Results/runs/<run-id>/` root.

Environment discovery, creation, verification, and failures are recorded under `.workflow/environment-bootstrap/` and summarized in `Scripts/Environment/Report.md`. Never report a profile as ready from package installation alone; its profile verification command must pass.

## Run log format

The run log is a canonical re-render, never a surgical insert. `build_report_text()` in `tools/_common.py` is the single renderer, so rebuilding the log is identical to appending to it by construction rather than by keeping two code paths in agreement. This also means a malformed file converges on the next write instead of accumulating damage.

```
<!-- SCMO-RUNLOG:START -->          one region, always last in the file
  <two-line bilingual prologue: tool-owned, do not edit inside>
  <!-- SCMO-RUN:<key>:START -->    one pair per run, used to dedupe by run
  ### `<run_id>` — complete
  - Status / 状态：`complete`
  - Checked / 检查时间：`2026-09-19T08:12:44+00:00`
  - Tasks / 任务：12/21 complete, 1 failed, 8 unfinished
  - Input signature / 输入签名：`<64-hex>`
  - Evidence / 证据：`.workflow/runs/run_001/run_summary.json`
  <!-- SCMO-RUN:<key>:END -->
<!-- SCMO-RUNLOG:END -->
```

- **Ownership.** Only the region between the `RUNLOG` markers is tool-owned. The `## Run log / 运行记录` heading, the note, and everything above the region is human text and is preserved verbatim. Human notes go above the region, never inside it.
- **Marker keys are not run ids.** A key must match `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; anything else (for example `a-->b`, which `plan_workflow.py` accepts) becomes `h` plus the first 16 hex digits of its SHA-256, because a key containing `-->` would terminate the HTML comment and destroy marker parsing permanently. The real run id still appears verbatim in the record body.
- **Parse order.** The region is found by locating `RUNLOG:START` and then searching for `RUNLOG:END` only *after* it, so an END before its START cannot happen. A run section missing its END is closed at the next START or at the region's end and re-emitted with both markers.
- **Counts come from evidence.** `complete` counts tasks whose `evidence_valid` is true, `failed` counts terminal scheduler failures, and the remainder is `unfinished` — a task that was never submitted and a task whose scheduler state says complete but whose signatures do not verify are both unfinished, not complete.
- **Union, disk wins.** Records already in the file survive a pruned or archived run directory, and the freshly written `run_summary.json` wins on collision. The log therefore accumulates. Ordering is by parsed `checked_at`, oldest first, with unparseable stamps last in file order; re-inspecting a run re-stamps it and moves it to the end.
- **The format is frozen at `run-log format 1`.** Change the number only for a change a reader must be warned about, and update this section with it.
