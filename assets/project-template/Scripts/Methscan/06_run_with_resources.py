#!/usr/bin/env python3
"""Run one workflow stage and record portable child-process resource usage."""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    return args


def exit_code(wait_status):
    if os.WIFEXITED(wait_status):
        return os.WEXITSTATUS(wait_status)
    if os.WIFSIGNALED(wait_status):
        return 128 + os.WTERMSIG(wait_status)
    return 1


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(timezone.utc).astimezone().isoformat()
    started = time.monotonic()
    try:
        process = subprocess.Popen(args.command)
    except OSError as error:
        record = {
            "command": args.command,
            "started_at": started_at,
            "finished_at": datetime.now(timezone.utc).astimezone().isoformat(),
            "elapsed_seconds": time.monotonic() - started,
            "return_code": 127,
            "launch_error": str(error),
        }
        with args.output.open("w") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return 127
    _, wait_status, usage = os.wait4(process.pid, 0)
    return_code = exit_code(wait_status)
    process.returncode = return_code
    record = {
        "command": args.command,
        "started_at": started_at,
        "finished_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "elapsed_seconds": time.monotonic() - started,
        "return_code": return_code,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "max_rss_kb": usage.ru_maxrss,
        "input_blocks": usage.ru_inblock,
        "output_blocks": usage.ru_oublock,
        "voluntary_context_switches": usage.ru_nvcsw,
        "involuntary_context_switches": usage.ru_nivcsw,
    }
    with args.output.open("w") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return return_code


if __name__ == "__main__":
    sys.exit(main())
