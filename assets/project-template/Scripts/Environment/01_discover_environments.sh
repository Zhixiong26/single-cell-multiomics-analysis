#!/usr/bin/env bash
# Print a read-only TSV inventory of active tools and known Conda environments.
# Set SCMO_ENV_ROOTS (colon-separated) to point at additional environment roots.
set -euo pipefail

printf 'record_type\tname\tpath\tversion\tavailable_tools\n'

tools=(python3 conda mamba samtools bismark bismark_methylation_extractor bowtie2 fastqc multiqc methscan allcools bgzip tabix)
for tool in "${tools[@]}"; do
  path=$(command -v "$tool" 2>/dev/null || true)
  [[ -n "$path" ]] || { printf 'active_tool\t%s\t\tmissing\t\n' "$tool"; continue; }
  version=$("$path" --version 2>&1 | sed -n '1p' || true)
  printf 'active_tool\t%s\t%s\t%s\t\n' "$tool" "$path" "${version//$'\t'/ }"
done
# $HOME/miniconda3/envs and $HOME/miniforge3/envs are conventional install
# locations, not assumptions about this site: they carry no project identity and
# expand per user. SCMO_ENV_ROOTS and CONDA_ENVS_PATH are searched in addition,
# which is how a project points at prefixes elsewhere on the filesystem.
roots=("$HOME/miniconda3/envs" "$HOME/miniforge3/envs")
if [[ -n "${SCMO_ENV_ROOTS:-}" ]]; then IFS=: read -r -a configured <<< "$SCMO_ENV_ROOTS"; roots+=("${configured[@]}"); fi
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
