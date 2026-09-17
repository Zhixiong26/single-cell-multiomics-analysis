# Intake and configuration

`init_project.py` accepts a mapping with `schema_version: 1`, `project_id`, optional organism/mitochondrial prefixes, references, samples, environments, scheduler settings, annotation, and analysis overrides. JSON-compatible YAML works without PyYAML.

Each sample supplies exactly: `sample_id`, `condition`, `batch`, `include`, `rna_path`, `rna_format`, `allc_root`, `allc_glob`, and `cell_id_prefix`. Biological metadata must be explicit. `rna_format` is blank, `10x_mtx`, `10x_zip`, or `10x_h5`. An included sample must provide RNA and/or ALLC.

ALLC roots default to recursive `**/*.allc.tsv.gz`; ALLCools-style `*_allc.gz` is also supported through an explicit `allc_glob`. Every selected compressed ALLC requires its adjacent `.tbi`. Use a narrower glob where unrelated files share the directory.

Generated interfaces are `config/project.yaml`, `samples.tsv`, `analysis.yaml`, `environments.tsv`, and `scheduler.yaml`. Relative data paths resolve against the generated project root. References may include `genome_id`, `chrom_sizes`, `blacklist`, optional `tss_bed`, and corresponding `_sha256` fields. The skill supplies no implicit genome or blacklist.

RNA enables Scanpy. ALLC enables VMR, ALLCools, and related MethylVI routes. Cell-type DMR and VMR+DMR require a readable annotation table. Missing modalities disable their routes rather than invalidating supported routes.
