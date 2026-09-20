# Documentation and evidence

Generated projects and major stages have bilingual README and Report files. README holds stable purpose, inputs, outputs, environments, parameters, order, entry points, and completion criteria. Report holds observed commands, jobs, resources, signatures, QC, failures, recovery, limitations, and unresolved decisions.

Two levels exist. The root `README.md`/`Report.md` are the project's contract and run log: the README is written once by `init_project.py` from the intake and the resolved configuration, and the Report accumulates one concise record per run. The module documents under `Scripts/<Module>/` (and `Scripts/Methylvi/<route>/`) carry the same two kinds of content scoped to one stage, and each is owned separately:

- The **stage README and the prose of the stage Report** are the stage's operating contract, ported from the template and bound to this project once by `init_project.py`. Their worked examples and numbers describe the reference project the template came from, and are labelled as such by the generated context block; they are never rewritten to claim they describe this project's data. A stage's stable purpose, inputs, outputs, and completion criteria are the same in every project, so a per-project copy of them is a copy, not a state record.
- The **stage Report's run-log region** is machine-owned and current. `refresh_stage_run_logs()` in `tools/_common.py` re-renders it from the same `.workflow/runs/*/run_summary.json` evidence that fills the root log, filtered to the tasks that belong to that stage, and `inspect_run.py` refreshes it on every inspection just as it refreshes the root log.

A run therefore keeps the root Report and each stage Report it touched in agreement by construction: both are rendered from one set of summaries by one renderer, and a stage whose tasks were not part of a run gets no record for it.

The root README carries one section the module documents do not: the `Data/` entries and where each came from. It is generated with the project and frozen with the rest of the file, so it describes the inputs as configured rather than as they were later moved. The root Report repeats that table above the run-log region together with two things nothing else states — which bundled reference files generation filled in, and any sample whose input was configured as a path outside the project. Both are above the region and are therefore human text the tools never rewrite.

Keep dynamic state under `.workflow/runs/<run-id>/`. A stage is complete only when summaries, artifacts, and completion markers agree; Slurm `COMPLETED` alone is insufficient. Figures require adjacent tabular or JSON QC/manifests.

Never rewrite the root README with runtime state: it is frozen after generation, because a document that reports progress stops being trustworthy as a description of what the project is. The same rule governs the stage READMEs and the prose above every stage run log. Runtime facts belong in the run-log regions. `update_report.py` rebuilds the root log and every stage log, and `inspect_run.py` refreshes the inspected run's record in both as a side effect of computing the run's verdict.

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
- **Stage logs are the same format, filtered.** A stage Report carries its own region, written by the same `build_report_text()` with a `stage` argument: records keep the identical shape and the identical marker keys, but only the tasks mapped to that stage by `stage_for_task()` appear in the counts and the section is omitted entirely for a run that never touched the stage. Because the key is the run id either way, re-inspecting a run updates its record in the root log and in each stage log it belongs to, and the two can never disagree about that run.
- **The format is frozen at `run-log format 1`.** Change the number only for a change a reader must be warned about, and update this section with it.
