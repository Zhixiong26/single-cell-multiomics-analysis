# Template maintenance

Keep algorithm scripts in stage directories, remove site paths and sample labels, and route project values through configuration or explicit CLI/environment inputs. Record source-to-template mapping in `migration-inventory.tsv`. Replace static Slurm resources with scheduler profiles.

Do not package logs, results, checkpoints, notebook output, credentials, Python caches, legacy uppercase `Config`, or dynamic state. The project generator must ignore `__pycache__`, bytecode, test caches, notebook checkpoints, OS metadata, and stale pre-schema-v2 configuration directories even when they exist in the installed skill checkout. Increment `VERSION` for releases and keep schema compatibility within a minor version.

Only runtime-reachable scripts belong in the execution signature. Interactive notebooks are templates, not executable workflow inputs, and are excluded. Do not retain independent Slurm wrappers that bypass the run DAG and its completion evidence.
