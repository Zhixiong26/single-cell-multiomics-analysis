#!/usr/bin/env bash
set -euo pipefail
printf 'record_type\tname\tpath\tversion\tavailable_tools\n'
tools=(python3 conda mamba samtools methscan allcools bgzip tabix)
for tool in "${tools[@]}"; do
  path=$(command -v "$tool" 2>/dev/null || true)
  [[ -n "$path" ]] || { printf 'active_tool\t%s\t\tmissing\t\n' "$tool"; continue; }
  version=$("$path" --version 2>&1 | sed -n '1p' || true)
  printf 'active_tool\t%s\t%s\t%s\t\n' "$tool" "$path" "${version//$'\t'/ }"
done
roots=("$HOME/miniconda3/envs" "$HOME/miniforge3/envs")
if [[ -n "${CONDA_ENVS_PATH:-}" ]]; then IFS=: read -r -a extra <<< "$CONDA_ENVS_PATH"; roots+=("${extra[@]}"); fi
for root in "${roots[@]}"; do
  [[ -d "$root" ]] || continue
  for env in "$root"/*; do
    [[ -d "$env" ]] || continue
    found=""
    for tool in "${tools[@]}"; do [[ -x "$env/bin/$tool" ]] && found="${found}${found:+,}${tool}"; done
    version=$([[ -x "$env/bin/python" ]] && "$env/bin/python" --version 2>&1 | sed -n '1p' || echo none)
    printf 'conda_env\t%s\t%s\t%s\t%s\n' "${env##*/}" "$env" "$version" "$found"
  done
done
