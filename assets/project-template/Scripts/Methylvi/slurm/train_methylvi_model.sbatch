#!/usr/bin/env bash
# Train and plot one of the configured method/version/feature models.
# Usage: sbatch [resource options] train_methylvi_model.sbatch allcools none 10000
#        sbatch [resource options] train_methylvi_model.sbatch vmr 0.01 30000
# Resources and log paths are injected by the submitter, not declared here.
#SBATCH --job-name=scmo_mvi_train

set -Eeuo pipefail
if (( $# != 3 )); then echo "Usage: sbatch $0 allcools|vmr VERSION FEATURES" >&2; exit 2; fi
method=$1
version=$2
features=$3
case "$method" in allcools|vmr) ;; *) echo "Unsupported method: $method" >&2; exit 3;; esac
script_dir=${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}/Scripts/Methylvi
source "$script_dir/00_experiment_config.sh"
[[ " $SCMO_FEATURE_TARGETS " == *" $features "* ]] || { echo "Unsupported feature target: $features" >&2; exit 4; }
if [[ "$method" == vmr ]]; then
  [[ " $SCMO_VMR_THRESHOLDS " == *" $version "* ]] || { echo "Unsupported VMR threshold: $version" >&2; exit 5; }
fi
mvi_root="$SCMO_MVI_EXPERIMENT_ROOT"
if [[ "$method" == allcools ]]; then
  route="$mvi_root/allcools/features_${features}"
else
  route="$mvi_root/vmr/var_${version}/features_${features}"
fi
input="$route/input.h5mu"
results="$route/results"
[[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 6; }
export SCMO_MVI_INPUT="$input"
export SCMO_MVI_RESULTS="$results"
# The VMR plotting entry point reads its default output location from this
# route-specific variable; keep it aligned with the shared MethylVI results.
export VMR_MVI_RESULTS="$results"
export SCMO_BATCH_KEY="${SCMO_BATCH_KEY:-sample_id}"
export SCMO_CELLTYPE_KEY="${SCMO_CELLTYPE_KEY:-cell_type}"
export SCMO_THREADS="${SLURM_CPUS_PER_TASK:-${SCMO_THREADS:-16}}"
export SCMO_SEED="${SCMO_SEED:-0}"
export SCMO_EPOCHS="${SCMO_EPOCHS:-500}"
export SCMO_BATCH_SIZE="${SCMO_BATCH_SIZE:-32}"
export SCMO_SUPERVISED_TARGET_KEY="${SCMO_SUPERVISED_TARGET_KEY:-$SCMO_CELLTYPE_KEY}"
export SCMO_SUPERVISED_TARGET_WEIGHTS="${SCMO_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export SCMO_SUPERVISED_NEIGHBORS="${SCMO_SUPERVISED_NEIGHBORS:-15}"
export SCMO_SUPERVISED_MIN_DIST="${SCMO_SUPERVISED_MIN_DIST:-0.5}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "${SCMO_PROJECT_ROOT:?set SCMO_PROJECT_ROOT}"
python_exe="$SCMO_METHYLVI_ENV/bin/python"
"$python_exe" - "$input" "$features" "${SCMO_EXPECTED_CELLS:?set SCMO_EXPECTED_CELLS}" <<'PY'
import os, sys, mudata, numpy as np
x = mudata.read_h5mu(sys.argv[1], backed="r")["mCG"]
expected = int(sys.argv[2])
cells = int(sys.argv[3])
if x.n_vars != expected or x.n_obs != cells:
    raise SystemExit(f"Input shape mismatch: {x.shape}, expected ({cells}, {expected})")
required = {os.environ["SCMO_BATCH_KEY"], "condition", os.environ["SCMO_CELLTYPE_KEY"]}
missing = required - set(x.obs.columns)
if missing:
    raise SystemExit(f"Missing metadata: {sorted(missing)}")
if set(np.asarray(x.var["selection_rank"], dtype=int)) != set(range(1, expected + 1)):
    raise SystemExit("Invalid selection_rank")
print(f"Input preflight passed: {x.n_obs} cells x {x.n_vars} features")
PY
"$python_exe" "$script_dir/shared/04_train_methylvi.py"
if [[ "$method" == allcools ]]; then
  "$python_exe" "$script_dir/allcools/07_plot_methylvi_umap.py" --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
else
  "$python_exe" "$script_dir/vmr/03_plot_vmr_methylvi_umap.py" --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
fi
"$python_exe" "$script_dir/shared/05_plot_supervised_umap.py"
"$python_exe" "$script_dir/shared/07_plot_methylation_qc.py" --input "$input" --output-dir "$results/qc"
touch "$route/model.COMPLETE"
