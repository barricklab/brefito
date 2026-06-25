# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What brefito is

`brefito` is a Python CLI that wraps Snakemake workflows for microbial genome assembly, polishing, annotation, and mutation prediction. Each workflow is a standalone `.smk` file; `brefito` resolves sample data and invokes `snakemake --use-conda` with the appropriate file and config arguments.

## Installation

```bash
# Create the base conda environment (provides snakemake, snakedeploy, biopython)
conda env create -f environment.yml
conda activate brefito

# Install the package in editable mode
pip install -e .
```

Each workflow rule uses its own per-tool conda environment defined in `src/brefito/workflow/envs/`.

## Running workflows

```bash
# Run a workflow on all samples
brefito <workflow-name>

# Run on specific samples
brefito <workflow-name> sample1 sample2

# Dry-run (shows snakemake plan)
brefito <workflow-name> --dry-run

# Pass config options to snakemake
brefito <workflow-name> --config key=value

# Specify an alternate reference directory (overrides 'references' default)
brefito predict-mutations-breseq-myref   # sets references=myref
brefito predict-mutations-breseq -r myref

# Recover from a locked snakemake directory
brefito <workflow-name> --pick-lock

# Force-reinstall a conda environment (e.g. to pick up a new prerelease build of a pinned tool)
brefito <workflow-name> --reinstall breseq
```

## Architecture

### CLI layer (`src/brefito/brefito.py`)

`main()` parses arguments, resolves the workflow name (including the `-<reference>` suffix shorthand), then constructs and executes a `snakemake --use-conda` command pointing at the matching `.smk` file. For polish workflows, it also handles renaming the output assembly through versioned `.N.breseq/.N.medaka/…` copies so iterations accumulate rather than overwrite.

### Sample data model (`load-sample-info.smk`)

Every workflow begins with:
```python
try: sample_info
except NameError:
    include: "load-sample-info.smk"
```

`load-sample-info.smk` defines the `SampleInfo` class and instantiates `sample_info`. It loads from `data.csv` if present, otherwise auto-scans these directories:
- `nanopore-reads/` — `*.fastq.gz`
- `illumina-reads/` — `*.R1.fastq.gz`, `*.R2.fastq.gz`, `*.fastq.gz` (SE)
- `assemblies/` — `*.fasta`, `*.fna`, `*.fa`
- `references/` — genome/annotation files

The `data.csv` format has columns `sample`, `type`, `setting`. Valid types: `nanopore`, `illumina-PE` (paired, use `{1|2}` placeholder), `illumina-R1`, `illumina-R2`, `illumina-SE`, `reference`, or arbitrary option keys (e.g., `breseq_option`). SRA accessions are supported via `sra://SRRxxxxxx`.

`SampleInfo` methods used in rules:
- `get_sample_list()` — all sample names
- `get_nanopore_read_list(sample)` / `get_illumina_PE_read_base_list(sample)` — file lists
- `get_reference_list(sample)` / `get_reference_arguments(sample, '-r ')` — reference files
- `get_reference_prefix()` — the active reference directory name

### Workflow rules (`src/brefito/workflow/rules/`)

Each `.smk` file is a self-contained workflow. Key workflows:

| Category | Workflows |
|----------|-----------|
| Assembly | `autocycler-assemble`, `assemble-flye`, `assemble-unicycler` |
| Polishing | `polish-breseq`, `polish-medaka`, `polish-polypolish`, `polish-polca` |
| Mutation prediction | `predict-mutations-breseq`, `predict-mutations-minimap2-breseq` |
| Read QC | `trim-illumina-reads`, `trim-nanopore-reads`, `filter-nanopore-reads`, `evaluate-nanopore-reads` |
| Annotation | `annotate-genomes` |
| Comparison | `compare-syri`, `compare-assemblies-breseq`, `compare-mutations-breseq`, `compare-genome-diffs-breseq` |
| Reads/data | `download-data`, `data-to-sra`, `merge-reads`, `align-reads` |
| Analysis | `classify-kraken2`, `search-blast`, `predict-cnv-breseq`, `tabulate-ssrs-breseq`, `evaluate-coverage-breseq` |

`autocycler-assemble` runs 6 assemblers (canu, flye, miniasm, necat, metamdbg, raven) × 4 subsampled read sets, then uses autocycler cluster/trim/resolve/combine to produce a consensus assembly in `assemblies/<sample>.fasta`. It requires `--config genome_size=<bp>`.

### Helper CLI tools

Installed alongside `brefito` via `pyproject.toml`:
- `normalize_assembly` — reindex/rename/sort contigs in a FASTA (useful after assembly)
- `fastq_folder_to_data_csv` — generate a `data.csv` from a folder of FASTQ files
- `canu_trim` — trim canu assembly output
- `extract_mobile_elements` — extract IS elements from ISEScan output

### Config keys passed through to workflows

Config is passed as `--config key=value` and accessed in rules as `brefito_config` (uppercased). Common keys:
- `data_csv` — path to data CSV (default: `data.csv`)
- `references` — reference directory name (default: `references`)
- `samples` — underscore-comma-separated sample filter (set automatically from CLI positional args)
- `genome_size` — required by autocycler-assemble
- `BRESEQ_THREADS`, `BRESEQ_OPTIONS`, `NO_DEFAULT_BRESEQ_OPTIONS` — breseq tuning
