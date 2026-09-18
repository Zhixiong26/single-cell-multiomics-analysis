# Intake and configuration

`init_project.py` generates schema v2 projects and accepts v1 or v2 intake during migration. Generated loaders continue to read schema v1 with deprecation warnings. JSON-compatible YAML works without PyYAML.

Each schema-v2 sample supplies `sample_id`, `condition`, `batch`, `include`, RNA fields, `allc_root`, `allc_glob`, `cell_id_prefix`, and optional `allc_cell_id_regex`/`allc_cell_id_replacement`. Sample ownership always comes from its manifest row. Without a regex, common ALLC suffixes are removed and the prefix is prepended. An included sample must provide RNA and/or ALLC.

ALLC roots default to recursive `**/*.allc.tsv.gz`; ALLCools-style `*_allc.gz` is also supported through an explicit `allc_glob`. Every selected compressed ALLC requires its adjacent `.tbi`. Use a narrower glob where unrelated files share the directory.

Generated interfaces are `config/project.yaml`, `samples.tsv`, `analysis.yaml`, `environments.tsv`, and `scheduler.yaml`. Relative paths resolve against the project root. Methylation routes require `chrom_sizes` and blacklist; TSS is optional. Actual SHA-256 values are always recorded, and a declared `_sha256` must match.

Annotation uses `table`, `profile`, and `review_status`. `table` is a cell TSV; `profile` is a guarded Scanpy cluster-mapping YAML. Legacy `path` is interpreted only as `table` and emits a warning. Cell-type DMR requires `review_status: approved` and rejects placeholder labels.

The orchestrator Python must be 3.9 or newer. Stage-specific interpreters and executables are declared with absolute paths in `environments.tsv`. Known compatible environments are reused read-only. Missing profiles are created from the bundled specs by `bootstrap_environments.py`, verified, and then recorded with absolute paths; see [environment provisioning](environments.md).

`scheduler.profiles` defines floor/target/ceiling resources for every task profile. `limited_profiles` selects memory-intensive profiles constrained by `max_parallel`; `validation_workers` controls full ALLC validation. Scanpy defaults explicitly include HVG flavor, batch key, and Scrublet expected doublet rate in generated `analysis.yaml`.

RNA enables Scanpy. ALLC enables VMR, ALLCools, and unsupervised MethylVI; the workflow generates `Unassigned` metadata when annotation is absent. Cell-type DMR and VMR+DMR require approved annotation. Missing modalities disable only unsupported routes.
