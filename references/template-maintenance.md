# Template maintenance

Keep algorithm scripts in stage directories, remove site paths and sample labels, and route project values through configuration or explicit CLI/environment inputs. Record source-to-template mapping in `migration-inventory.tsv`. Replace static Slurm resources with scheduler profiles.

Do not package logs, results, checkpoints, notebook output, credentials, or dynamic state. Increment `VERSION` for releases and keep schema compatibility within a minor version.
