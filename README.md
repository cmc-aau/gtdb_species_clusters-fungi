# gtdb_species_clusters-fungi

> **Disclaimer:** This repository is intended for internal use by the GTDB development team. It is tailored to GTDB's internal data layout and release process, and it is not generally recommended for external use. Interfaces, input formats and outputs may change between releases without notice.

`gsc_fungi` is a toolkit for establishing, updating and validating GTDB fungal species clusters. It is a fungal adaptation of [GTDB Species Clusters](https://github.com/Ecogenomics/gtdb-species-clusters).

## Installation

Create a conda environment from one of the environment YAML files in the root of the repository. Each YAML file pins exact versions of all dependencies (including [skani](https://github.com/bluenote-1577/skani) and [gtdblib](https://pypi.org/project/gtdblib/)). For example, to use `gtdb_r237_fungi-environment.yml`:

```bash
mamba env create -f gtdb_r237_fungi-environment.yml
conda activate gtdb_r237_fungi
```

The toolkit isn't installed into the environment; it is run from source. With the environment active, put the repository on your `PATH` and `PYTHONPATH`:

| Variable | Add | Why |
|---|---|---|
| `PATH` | `<repo>/bin` | makes the `gsc_fungi` command available |
| `PYTHONPATH` | `<repo>/src` | lets Python import the `gsc_fungi` package |

where `<repo>` is the root of this repository. For example, from the root of the repository:

```bash
export PATH=$PWD/bin:$PATH
export PYTHONPATH=$PWD/src${PYTHONPATH:+:$PYTHONPATH}
gsc_fungi -h
```

Nothing else needs to be added. The environment provides all dependencies, including gtdblib and skani. Don't add a local clone of gtdb-lib to `PYTHONPATH` unless you need unreleased changes, because it would override the version pinned in the environment.

## Input parameters file

Most commands take a single positional argument: a TOML file listing the input files and output directories. [example_files/gtdb_inputs.toml](example_files/gtdb_inputs.toml) is an example of the expected format. It has these parts:

| Section | Contents |
|---|---|
| `release`, `root_dir` | GTDB release name, and the directory that all relative paths in the file are resolved against |
| `[data_files]` | genome paths, NCBI assembly summary, NCBI taxonomy, genome metadata and marker set statistics |
| `[nomenclature_files]` | MycoBank data and the lists of sanctioned names |
| `[ledgers]` | QC exception list and hybrid species list |
| `[output_dirs]` | output directory for each command, keyed by the command's name |

Each command writes its results and a `gsc_fungi.log` file to its entry in `[output_dirs]`. Later commands read the outputs of earlier commands from these directories, so the commands must be run in order.

## Commands

Run `gsc_fungi <command> -h` for the full help of any command. All commands accept `--silent` to suppress output.

### Create de novo species clusters

#### `cdn_qc_genomes`

Quality checks all genomes. A genome passes if it meets the marker gene thresholds below and isn't removed by one of the optional NCBI or hybrid filters. Genomes on the QC exception ledger are retained or removed regardless of these criteria.

Writes `qc_passed.tsv`, `qc_failed.tsv` and `species_fail_qc.tsv` (NCBI species that lose all their genomes during QC).

```bash
gsc_fungi cdn_qc_genomes gtdb_inputs.toml
```

| Parameter | Description |
|---|---|
| `input_params_file` | TOML file with input parameters |
| `--min_perc_markers` | minimum percentage of phylogenetic marker genes, 0 to 100 (default: 30) |
| `--max_perc_dup_markers` | maximum percentage of duplicated, lineage-specific, single-copy genes, 0 to 100 (default: 10) |
| `--remove_ncbi_exclusions` | remove genomes that NCBI excluded from RefSeq for selected reasons (default: off) |
| `--remove_hybrids_on_list` | remove genomes from species on the hybrid species ledger (default: off) |
| `--keep_ncbi_hybrid_species` | keep genomes from NCBI species with hybrid names such as `A x B` (removed by default) |

#### `cdn_select_sp_reps`

Selects a representative genome for each NCBI species, using genomes that passed QC. Genomes from the type strain of the species are preferred, followed by genomes from the type strain of a subspecies, and then NCBI representative genomes. If none of these exist, the genome with the highest assembly score is chosen. When a species has several type strain genomes and none is the NCBI representative, the representative is chosen from the genomes within 99% ANI of most of the type strain genomes.

Reads `qc_passed.tsv` from `cdn_qc_genomes`. Writes `sp_reps.tsv`, plus `manual_reps.tsv` listing species whose representative must be chosen by hand.

```bash
gsc_fungi cdn_select_sp_reps gtdb_inputs.toml
```

| Parameter | Description |
|---|---|
| `input_params_file` | TOML file with input parameters |

#### `cdn_cluster_named_sp`

Creates ANI-based species clusters around the representatives of named species. Pairwise ANI and alignment fraction (AF) are calculated with skani. Representatives within 95% ANI and 50% AF of each other are merged, and the name that is retained follows MycoBank naming priority, with sanctioned names always taking priority. The remaining genomes that passed QC are then assigned to the clusters.

Reads `qc_passed.tsv` from `cdn_qc_genomes` and `sp_reps.tsv` from `cdn_select_sp_reps`. Writes `sp_clusters-ani95.tsv`, along with tables describing which species were merged and how naming priority was resolved.

```bash
gsc_fungi cdn_cluster_named_sp gtdb_inputs.toml --cpus 32
```

| Parameter | Description |
|---|---|
| `input_params_file` | TOML file with input parameters |
| `-c`, `--cpus` | number of CPUs (default: 1) |

#### `cdn_cluster_de_novo`

Listed in the help menu but not yet implemented.

### Update species clusters

#### `u_qc_genomes`

Quality checks genomes when updating species clusters from a previous release. It uses the same criteria and writes the same files as `cdn_qc_genomes`, but the NCBI exclusion and hybrid filters are on by default. The TOML file needs a `u_qc_genomes` entry in `[output_dirs]`.

```bash
gsc_fungi u_qc_genomes gtdb_inputs.toml
```

| Parameter | Description |
|---|---|
| `input_params_file` | TOML file with input parameters |
| `--min_perc_markers` | minimum percentage of phylogenetic marker genes, 0 to 100 (default: 30) |
| `--max_perc_dup_markers` | maximum percentage of duplicated, lineage-specific, single-copy genes, 0 to 100 (default: 10) |
| `--skip_ncbi_exclusion` | don't remove genomes based on NCBI exclusion metadata |
| `--keep_hybrids` | don't remove hybrid genomes, whether from the hybrid species ledger or with NCBI hybrid names |

### Utilities

These commands don't use the TOML file. Their inputs are given directly on the command line.

#### `taxa_by_rank`

Counts the genomes and named taxa at each rank of the NCBI taxonomy, and the genomes (including MAGs and SAGs) missing a classification at each rank. Writes `taxa_by_rank.tsv` and `missing_classifications_by_rank.tsv`.

```bash
gsc_fungi taxa_by_rank \
    -i ncbi_r232_standardized.tsv \
    -a assembly_summary_fungi_genbank.txt \
    -o taxa_by_rank_out
```

| Parameter | Description |
|---|---|
| `-i`, `--ncbi_taxonomy_file` | standardized NCBI taxonomy file with subranks (required) |
| `-a`, `--assembly_summary_file` | GenBank assembly summary file from NCBI (required) |
| `-o`, `--out_dir` | output directory (required) |

#### `plot_genome_stats`

Plots histograms of common genome statistics: BUSCO completeness, contamination and quality, contig count, contig N50 and ambiguous bases. Writes each plot as PNG and SVG, plus the data in `qc_table.tsv`.

```bash
gsc_fungi plot_genome_stats \
    -b busco_summary.tsv \
    -a assembly_summary_fungi_genbank.txt \
    -o genome_stats_plots
```

| Parameter | Description |
|---|---|
| `-b`, `--busco_file` | file with BUSCO genome quality estimates (required) |
| `-a`, `--assembly_summary_file` | GenBank assembly summary file from NCBI (required) |
| `-o`, `--out_dir` | output directory (required) |

## License

GPL-3.0. See [LICENSE](LICENSE).
