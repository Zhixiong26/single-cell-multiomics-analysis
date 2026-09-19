# Intake and configuration

`init_project.py` generates schema v2 projects and accepts v1 or v2 intake during migration. Generated loaders continue to read schema v1 with deprecation warnings. JSON-compatible YAML works without PyYAML.

Each schema-v2 sample supplies `sample_id`, `condition`, `batch`, `include`, RNA fields, `allc_root`, `allc_glob`, `cell_id_prefix`, and optional `allc_cell_id_regex`/`allc_cell_id_replacement`. Sample ownership always comes from its manifest row. Without a regex, common ALLC suffixes are removed and the prefix is prepended. An included sample must provide RNA and/or ALLC.

ALLC roots default to recursive `**/*.allc.tsv.gz`; ALLCools-style `*_allc.gz` is also supported through an explicit `allc_glob`. Every selected compressed ALLC requires its adjacent `.tbi`. Use a narrower glob where unrelated files share the directory.

Generated interfaces are `config/project.yaml`, `samples.tsv`, `analysis.yaml`, `environments.tsv`, and `scheduler.yaml`. Relative paths resolve against the project root. Methylation routes require `chrom_sizes` and blacklist; TSS is optional. Actual SHA-256 values are always recorded, and a declared `_sha256` must match.

## Input data

Analysis inputs are read only from the project's own `Data/` directory, so the project stays reproducible from its own tree. An optional `data` block declares where each input comes from, one entry per `Data/<name>`:

```yaml
data:
  sources:
    - {name: Matrix, path: /shared/lung-tissue/Matrix}
    - {name: ctrl_10x, url: "https://example.org/ctrl.tar.gz", sha256: "<64 hex>"}
```

`name` is a single path component. Exactly one of `path` (link to a location that already holds the data) or `url` (fetch it) is required. A relative `path` is anchored to the **intake file's directory** rather than the working directory, so one intake links the same targets wherever it is run. A `url` source must carry the sha256 of the finished file: every later stage reads whatever sits at the destination as data, and a truncated transfer is not distinguishable from a complete one by size.

`init_project.py` creates the links during generation and records each entry's resolved target in `config/project.yaml`, so the link can be rebuilt and the origin is part of the input signature. It does not download: a `url` source is left pending, and `tools/link_data.py --project PROJECT --execute` fetches it into a dot-prefixed sibling file, verifies the digest, and only then moves it into place. `link_data.py` without `--execute` reports the plan and changes nothing. It never overwrites an entry it did not create — an existing file or a link to somewhere else means someone has to decide which copy is real. `--only NAME` restricts a run to named entries.

A path that is not visible from the generating host still produces the link; the target may only be mounted in the execution context. `validate_project.py` refuses a link that does not resolve, naming it, so a mount visible on a login node and absent on a compute node fails where the data would have been read rather than inside a stage. `--allow-missing-paths` skips those checks for validating a configuration somewhere the data is not.

The root README lists every `Data/` entry with its origin, and the root Report records the same at generation, including any sample input that was configured as a project-external path instead of through `Data/`.

## Reference files

`references.chrom_sizes` and `references.blacklist` may be omitted: the skill ships `Supplementary/hg38.canonical.chrom.sizes` and `Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz`, copied into every generated project, and `init_project.py` fills in whichever of the two the intake leaves blank and `blacklist_md5` alongside the bundled blacklist. It records that it did so in the generated Report. The bundled files are human hg38, so they are applied only to a human project that declares a methylation route and either no `genome_id` or an hg38 one; a project declaring another organism or genome keeps the preflight failure and is told why, rather than silently inheriting coordinates for the wrong genome. When the intake declares no genome, generation also writes `references.genome_id: GRCh38`, so the recorded genome cannot contradict the files it selected; the Report says the key was inferred rather than declared, because the frozen README renders it as a project fact. Declaring either key in the intake overrides the default for that key alone, and a self-supplied blacklist needs its own `blacklist_md5`.

Annotation uses `table`, `profile`, and `review_status`. `table` is a cell TSV; `profile` is a guarded Scanpy cluster-mapping file (JSON or YAML). Legacy `path` is interpreted only as `table` and emits a warning. Cell-type DMR requires `review_status: approved` and rejects placeholder labels.

A profile binds a reviewed mapping to the run it was reviewed from, so it cannot be transplanted onto different parameters or a different cluster set. It holds `analysis_signature` (must equal the signature of the run being planned), `expected_clusters` (must equal that run's cluster set exactly), and `annotations`, which maps each cluster to either a bare label or `{cell_type, confidence, evidence}`. The Scanpy route writes the file for you: `tools/record_annotation_review.py --run-id ID --worksheet review.tsv` prepares a per-cluster worksheet seeded with each cluster's majority label, and the same command with `--mapping review.tsv` records the reviewed result, refusing placeholder labels and any cluster set that does not match the run. The cell-type DMR routes stay unplannable until `review_status` is `approved` as well.

The orchestrator Python must be 3.9 or newer. Stage-specific interpreters and executables are declared with absolute paths in `environments.tsv`. Known compatible environments are reused read-only. Missing profiles are created from the bundled specs by `bootstrap_environments.py`, verified, and then recorded with absolute paths; see [environment provisioning](environments.md).

`scheduler.profiles` defines floor/target/ceiling resources for every task profile. `limited_profiles` selects memory-intensive profiles constrained by `max_parallel`; `validation_workers` controls full ALLC validation. Scanpy defaults explicitly include HVG flavor, batch key, and Scrublet expected doublet rate in generated `analysis.yaml`.

`analysis.methscan.vmr_thresholds` must be a non-empty list of unique numeric values in `(0, 1]`. `analysis.methylvi.feature_targets` must be a non-empty list of unique positive integers. Invalid empty lists are rejected during preflight rather than failing inside DAG construction.

RNA enables Scanpy. ALLC enables VMR, ALLCools, and unsupervised MethylVI; the workflow generates `Unassigned` metadata when annotation is absent. Cell-type DMR and VMR+DMR require approved annotation. Missing modalities disable only unsupported routes. Explicit route subsets validate only the input modalities, references, and environments used by the resulting closed DAG; `auto` validates all automatically selected routes.

## Worked example

`config/examples/` ships a completed two-donor intake (`example-ipf-tissue.yaml`), the same rows in generated-manifest form (`example-samples.tsv`), and a README recording that study's observed magnitudes and the values that must be re-derived rather than copied. Generated projects never read this directory; it exists so a new project has one filled-in instance of every field. The scheduler block there is deliberately empty — partitions, accounts, and node lists describe your cluster, not someone else's.

## Variable namespace

Stage scripts and wrappers read project values from the neutral `SCMO_*` namespace, which the run DAG injects from `config/`. Route configs may keep a source-flavoured inner name (`MVI_*`, `VMR_*`, bare lowercase) internally, but each one re-exports the resolved value under its `SCMO_*` name, and new code reads only the neutral name. Exported `SCMO_*` values are never silently overwritten by a later-sourced config; alias assignments use `:-` guards.

Paths in the packaged scripts are project-root-relative or environment-supplied. No project root, login or compute node name, shared filesystem path, scheduler account, or absolute interpreter path is hardcoded outside `config/examples/`.
