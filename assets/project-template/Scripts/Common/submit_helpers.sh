#!/usr/bin/env bash
# Sourced by the standalone submit_*.sh wrappers.
#
# Stage wrappers carry no site-specific #SBATCH directives. scmo_submit reads the
# resource values for a named profile from config/scheduler.yaml and passes them,
# together with the log paths and the dependency, on the sbatch command line.
# This is the same contract tools/submit_workflow.py uses for the managed DAG.

scmo_submit() { # PROFILE DEPENDENCY JOBFILE [ARG ...]
  local profile=$1 dependency=$2 jobfile=$3; shift 3
  local root=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}
  local name flags log_dir
  name=$(basename "$jobfile" .sbatch)
  log_dir=${SCMO_SUBMIT_LOG_DIR:-$(dirname "$jobfile")/logs}
  mkdir -p "$log_dir"
  flags=$(bash "$root/Scripts/Common/scheduler_flags.sh" "$profile")
  local -a cmd=(sbatch --parsable --job-name "scmo_${name}"
    --output "$log_dir/${name}_%j.out" --error "$log_dir/${name}_%j.err")
  local -a resource_flags=()
  if [[ -n "$flags" ]]; then read -r -a resource_flags <<< "$flags"; fi
  cmd+=("${resource_flags[@]}")
  if [[ -n "$dependency" ]]; then cmd+=(--dependency "$dependency"); fi
  cmd+=("$jobfile" "$@")
  "${cmd[@]}"
}
