# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`gsc_fungi` is a toolkit for establishing, updating, and validating GTDB fungal species clusters. It is a fungal adaptation of GTDB Species Clusters (https://github.com/Ecogenomics/gtdb-species-clusters; local clone at `~/src/gtdb-species-clusters`). When asked to do something "the way GTDB Species Clusters does it", check that repo's `origin/master`, which may be ahead of the local checkout.

## Running

There is no packaging (`setup.py`/`pyproject.toml`), no test suite, and no linter configuration. The code runs from source:

```bash
# dependencies: gtdblib (~/src/gtdb-lib), rich, numpy, matplotlib, unidecode; skani on PATH for ANI steps
conda activate gtdb_species_clusters-fungi-r232     # or call ~/.conda/envs/gtdb_species_clusters-fungi-r232/bin/python
export PYTHONPATH=~/src/gtdb-lib:$PWD/src:$PYTHONPATH

bin/gsc_fungi <command> -h                          # wrapper for: python3 src/gsc_fungi "$@"
bin/gsc_fungi cdn_qc_genomes example_files/gtdb_inputs.toml
```

To check a change without writing into the real release directories, copy the TOML and point `output_dirs.<command>` (and any ledger being tested) at absolute paths in a scratch directory. `os.path.join` with an absolute path ignores `root_dir`, so inputs are still read from the real data. QC on release232 (~22k genomes) takes seconds.

## Architecture

**CLI flow.** `src/gsc_fungi/__main__.py` defines argparse subcommands, sets up the `rich` logger, adds a file handler writing `gsc_fungi.log` into the command's output directory, and logs the version and command line. It then calls `OptionsParser.run()` in `main.py`. Each handler there reads the TOML, calls `log_command_inputs()` to log the TOML entries and earlier-step outputs it uses, checks that the files exist, and then calls the `run()` method of a step class (`QcGenomes`, `SelectSpeciesRepresentative`, `ClusterNamedSpecies`, …). Use `logging.getLogger('rich')` so messages reach the log file.

**TOML input parameters.** Nearly all inputs come from one TOML file (`input_params_file`), with sections `[data_files]`, `[nomenclature_files]`, `[ledgers]` and `[output_dirs]`. Paths are relative to `root_dir` and are resolved lazily by `OptionsParser.args_abs_path()`, so `params` keeps the raw values. Each command's output directory is the `[output_dirs]` entry named after the subcommand. `example_files/gtdb_inputs.toml` is the current template. `docs/example_cdn_inputs.toml` is outdated.

**Pipeline.** `cdn_*` commands build clusters de novo, and each step reads files written by earlier steps from their output directories:
1. `cdn_qc_genomes` writes `qc_passed.tsv`, `qc_failed.tsv` and `species_fail_qc.tsv`.
2. `cdn_select_sp_reps` reads `qc_passed.tsv` and writes `sp_reps.tsv`.
3. `cdn_cluster_named_sp` reads both earlier outputs, runs skani ANI/AF, and applies MycoBank naming priority with sanctioned names.

`cdn_cluster_de_novo` is listed in the help menu but not implemented.

`u_qc_genomes` (update workflow) shares the `qc_genomes()` handler but keeps the older `--skip_ncbi_exclusion` and `--keep_hybrids` flags. These are translated in `main.py` to the same `QcCriteria`.

**QC (`qc_genomes.py`).** A genome passes if it meets all of these:
- `Fun77_SingleCopy_pct_AAI100 >= --min_perc_markers` (default 30)
- `SpecSet_Duplicated_pct_AAI99 <= --max_perc_dup_markers` (default 10)
- it isn't filtered by an optional NCBI or hybrid filter

The marker columns come from `marker_set_stats_file`. Genomes on the QC exception ledger always pass. The optional filters are:
- `--remove_ncbi_exclusions` removes genomes whose `excluded_from_refseq` includes a reason in `NCBI_EXCLUSION_FILTERING_CRITERIA` (`ncbi.py`). It is off by default.
- `--remove_hybrids_on_list` removes species on the hybrid ledger. It is off by default.
- Species whose NCBI name contains ` x ` (`is_hybrid_species()`) are removed by default. `--keep_ncbi_hybrid_species` keeps them.

The order of the per-genome `failed_*` flags written in `qc_genomes()` must match the headers of `qc_failed.tsv` and `species_fail_qc.tsv`. Keep them in sync when adding or removing a criterion.

**Ledgers.** Columns are found by header name, and files may be gzipped:
- QC exception list: `genome_id`, `exception_reason`, `route`
- hybrid species list: `species`, `hybrid_class`

Ledger formats change between releases, so check the headers of the actual files under `/projects/PC/GTDB/gtdb_fungal/metadata/<release>/ledgers/` rather than relying on `example_files/`.

**Gotchas.**
- `ClusterNamedSpecies.calculate_pairwise_ani()` has an `if True: # *** DEBUGGING` guard, so ANI is always recomputed and the `ani_af.pkl` cache is never read.
- `cluster_named_sp.run()` only loops over ANI threshold 95.
- Thresholds live in `defaults.py`.
- `scripts/` holds standalone one-off analysis scripts with hard-coded release paths. They are not part of the CLI.
- Release notes go in `src/gsc_fungi/VERSION`, newest first. The first line is the version string read by `__init__.py`.

## Pull requests

Every PR must bump the version in `src/gsc_fungi/VERSION`. To do this, add a new entry at the top of the file: the new version number on the first line, followed by a brief `- ` bullet for each change and a blank line before the previous entry. The first line of the file is what `__init__.py` reads as `__version__`.
