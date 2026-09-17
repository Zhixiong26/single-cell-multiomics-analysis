#!/usr/bin/env bash
set -Eeuo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python_exe=${PYTHON_EXE:-python3}
"$python_exe" -m unittest discover -s "$root/tests" -p 'test_*.py' -v
