# Environment discovery and provisioning

Environment handling has two phases. First inspect declared paths and discover compatible Conda prefixes without changing them. Then create only the unresolved profiles as isolated environments.

The generated project carries versioned specs for three profiles:

| Profile | Stages | Verification |
|---|---|---|
| `analysis_core` | orchestrator, Scanpy, ALLCools | imports Scanpy, ALLCools, Harmony, Scrublet, Leiden and YAML |
| `methscan` | MethSCAn VMR/DMR | runs `methscan --version` |
| `methylvi` | MethylVI training and plotting | imports `scvi.external.METHYLVI`, PyTorch, Scanpy and MuData |

Run a read-only plan first:

```bash
PYTHON=/path/to/python-3.9-or-newer
"$PYTHON" tools/bootstrap_environments.py --project PROJECT
```

If profiles remain unresolved, create them and update `config/environments.tsv`:

```bash
"$PYTHON" tools/bootstrap_environments.py \
  --project PROJECT \
  --execute
```

The default prefixes are `PROJECT/.environments/analysis-core`, `methscan`, and `methylvi`, with a writable package cache under the same root. Use `--prefix-root /new/writable/location` when project storage is unsuitable. `mamba` is preferred when present; otherwise `conda` is used. A specific manager may be passed with `--manager`.

Do not install into an active base environment, run `conda update`, or modify a discovered shared prefix. The bootstrapper refuses a non-empty target that fails verification. A failed creation is retained for diagnosis rather than silently deleted. Plans and results are saved under `.workflow/environment-bootstrap/`.

Environment creation requires package-channel/network access and sufficient storage. If no Conda-compatible manager is available, installation fails before analysis with a clear prerequisite instead of bypassing validation.
