#!/usr/bin/env python3
"""Update the generated bilingual Report from machine-readable run evidence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from _common import WorkflowError


START = "<!-- SCMO-RUN-STATE:START -->"
END = "<!-- SCMO-RUN-STATE:END -->"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    summary_path = root / ".workflow" / "runs" / args.run_id / "run_summary.json"
    if not summary_path.is_file():
        raise WorkflowError("run_summary.json is absent; run inspect_run.py first")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = [
        START,
        "## Latest run / 最近运行",
        "",
        "- Run / 运行：`%s`" % summary["run_id"],
        "- Status / 状态：`%s`" % summary["status"],
        "- Checked / 检查时间：`%s`" % summary["checked_at"],
        "",
        "| Task / 任务 | Job | State / 状态 |",
        "|---|---:|---|",
    ]
    for item in summary["tasks"]:
        rows.append("| %s | %s | %s |" % (item["task"], item["job_id"], item["state"]))
    rows.extend([END, ""])
    block = "\n".join(rows)
    report = root / "Report.md"
    text = report.read_text(encoding="utf-8") if report.is_file() else "# Analysis report / 分析报告\n\n"
    if START in text and END in text:
        text = text[:text.index(START)] + block + text[text.index(END) + len(END):]
    else:
        text = text.rstrip() + "\n\n" + block
    report.write_text(text, encoding="utf-8")
    print(str(report))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, OSError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
