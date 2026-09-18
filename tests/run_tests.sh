#!/usr/bin/env bash
set -Eeuo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python_exe=${PYTHON_EXE:-python3}
"$python_exe" -c 'import sys; assert sys.version_info >= (3, 9), "tests require Python >= 3.9; set PYTHON_EXE to the verified orchestrator Python"'
"$python_exe" -m unittest discover -s "$root/tests" -p 'test_*.py' -v
