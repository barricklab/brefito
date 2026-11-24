# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git: never run git commands without explicit approval

Do NOT run ANY git command that changes state — `git commit`, `git add`, `git merge`,
`git push`, `git rebase`, `git reset`, `git checkout`/`switch`, `git stash`, `git branch`,
tag creation, etc. — unless the user has explicitly told you to in the current request.
Wait for an explicit OK each time; approval for one git action does not carry over to the
next. Read-only git commands (`git status`, `git diff`, `git log`, `git show`) are fine to
run when helpful. When you believe a commit/merge/push is warranted, propose it and ask
first.

## Working directory: use the worktree, not the main checkout

When you are invoked in a git worktree (e.g. under `.claude/worktrees/<name>`), do ALL of
your work inside that worktree. Edit files by their path under the worktree directory, and
run git and other commands from there.

Do NOT modify the main checkout (or commit/merge onto its checked-out branch) unless the
user explicitly directs you to. In particular:

- Watch for absolute paths: a path under `.claude/worktrees/<name>/` points at the worktree,
  while the same repo-relative path outside that prefix points at the MAIN checkout. Always
  edit the worktree path so your changes land on the worktree branch.
- If a plan or prior message references a main-checkout path, translate it to the equivalent
  worktree path before editing.
- Commit on the worktree branch. Only touch `main` or the main checkout when the user asks
  you to merge or commit there.

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
- `nanopore-reads/`, `pacbio-clr-reads/`, `pacbio-hifi-reads/` — `*.fastq.gz`
- `illumina-reads/` — `*.R1.fastq.gz`, `*.R2.fastq.gz`, `*.fastq.gz` (SE)
- `assemblies/` — `*.fasta`, `*.fna`, `*.fa`
- `references/` — genome/annotation files

The `data.csv` format has columns `sample`, `type`, `setting`. Valid types: `nanopore`, `pacbio-clr` (PacBio RS II/Sequel continuous long reads), `pacbio-hifi` (PacBio HiFi/CCS), `illumina-PE` (paired, use `{1|2}` placeholder), `illumina-R1`, `illumina-R2`, `illumina-SE`, `reference`, or the name of any `--config` option (see below). SRA run accessions are supported via `sra://SRRxxxxxx` (or `ERR…`/`DRR…`, the ENA and DDBJ run accessions, which the SRA also serves).

**Long reads.** `nanopore`, `pacbio-clr` and `pacbio-hifi` are the long-read types (`LONG_READ_TYPES` in `load-sample-info.smk`). A sample may have reads of at most one of them; mixing them in one sample is an error at load time. Each type has its own directories, `<type>-reads/`, `<type>-reads-trimmed/`, `<type>-reads-filtered/`, `<type>-reads-subsampled/`, so nanopore paths are unchanged. Aligned reads are `nanopore_reads.*`, `pacbio_clr_reads.*` or `pacbio_hifi_reads.*`. Rules take the read type from an `{lrt}` wildcard (or `{long_read_prefix}` in align-reads), and per-type tool settings come from `LONG_READ_SETTINGS` via `long_read_setting(type, key)`:

| | nanopore | pacbio-clr | pacbio-hifi |
|---|---|---|---|
| minimap2 (align-reads) | `map-ont` | `map-pb` | `map-hifi` |
| minimap2 overlaps (miniasm) | `ava-ont` | `ava-pb` | `ava-pb` |
| flye | `--nano-raw` | `--pacbio-raw` | `--pacbio-hifi` |
| canu | `-nanopore` | `-pacbio` | `-pacbio-hifi` |
| autocycler helper `--read_type` | (default, `ont_r10`) | `pacbio_clr` | `pacbio_hifi` |
| filtlong default | `--min_length 1000 --keep_percent 90` | same | `--min_length 1000` |
| breseq default (minimap2-breseq) | `-x --require-match-fraction 0.0` | `--require-match-fraction 0.0` | `--require-match-fraction 0.0` |

ONT-only steps: porechop/porechop_ABI (`NANOPORE_TRIMMING` applies to nanopore reads only; PacBio reads are linked through untrimmed), medaka (`polish-medaka` stops with an error for PacBio samples), NECAT (left out of autocycler and trycycler for PacBio), and metaMDBG has no CLR mode (left out of autocycler for `pacbio-clr`). `NANOPORE_FILTERING` and `FILTLONG_OPTIONS` apply to every long-read type. Rule and workflow names keep their historical `nanopore` names (`filter-nanopore-reads`, `align_nanopore_reads`, ...) although they handle all long reads.

A `sample` containing `*`, `?` or `[...]` is a shell glob (no regex). Glob rows are expanded after the whole file is read, so they apply to every sample named exactly anywhere in the file, not just earlier rows. A glob matching nothing is warned about.

**Per-sample options in data.csv.** A row whose `type` is a config key (case-insensitive) sets that option for the matching sample(s), e.g. `REL1166A,BRESEQ_OPTIONS,--polymorphism-mode`. Resolution for key `K` (`get_config_value()` in `load-sample-info.smk`):
1. Workflow default.
2. data.csv `K`: a row naming the sample exactly beats a glob row; among rows at the same level, the later row wins. Replaces the default.
3. data.csv `K_ADD`: every matching row is appended, in file order.
4. Command-line `--config K=...` replaces all of the above, with a WARNING listing the affected samples.
5. Command-line `--config K_ADD=...` or `brefito --config-add K=...` is put in front of the result (command line first, then data.csv).

On/off flags (`config_is_true(name, default, sample=None)`): a blank data.csv setting means on; `0/false/no/off` means off. Keys brefito sets itself (`SAMPLES`, `REFERENCES`, `DATA_CSV`, `BREFITO_PACKAGE_PATH`, `QUERY_FILE_PATH`, `BRESEQ_PRERELEASE`) cannot be set in data.csv.

Some options cannot vary between samples. These are rule-order choices like `ILLUMINA_TRIMMING`/`NANOPORE_TRIMMING`/`NANOPORE_FILTERING`, thread counts, and keys that name output files (`*_REGION`, `*_FORMAT`, `SORT`/`MERGE`/`INDEX`). Per-read-file options like `FASTP_*` also cannot vary between the samples that share a read file. Workflows read these with `get_uniform_config_value()`, which accepts a data.csv value only if every relevant sample agrees and otherwise exits with an error.

When writing workflows, read per-sample options in a rule's `params:` (`lambda wildcards: get_config_value('KEY', DEFAULT, wildcards.sample)`) and use `{params.x}` in `shell:`, rather than a module-level global.

`SampleInfo` methods used in rules:
- `get_sample_list()` — all sample names
- `get_nanopore_read_list(sample)` / `get_illumina_PE_read_base_list(sample)` — file lists
- `get_long_read_type(sample)` (None, `nanopore`, `pacbio-clr` or `pacbio-hifi`), `get_long_read_list(sample)`, `get_long_read_base_list(sample)`, `get_samples_with_long_reads()`, `get_long_read_path(sample, stage, suffix)` — long reads of whichever type the sample has; use these, not the nanopore ones, in new workflows
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
| GenomeDiff files | `validate-genome-diff-gdtools`, `clean-genome-diffs-gdtools` |
| Reads/data | `download-data`, `data-to-sra`, `merge-reads`, `align-reads` |
| Analysis | `classify-kraken2`, `search-blast`, `predict-cnv-breseq`, `tabulate-ssrs-breseq`, `evaluate-coverage-breseq` |

`annotate-genomes` runs Prokka, then ISEScan, then PhiSpy, and combines them with
`breseq CONVERT-REFERENCE`. Each stage has a skip flag, settable globally or per sample in data.csv:
`SKIP_PROKKA` keeps the reference's own gene annotations (so the only additions are ISEScan's IS elements
and PhiSpy's prophages), `SKIP_ISESCAN` and `SKIP_PHISPY` drop those stages. The output path
`annotated-<references>/<sample>.gbk` is the same whichever are set. Adding IS element annotations to a
reference that lacks them is `--config SKIP_PROKKA=1 SKIP_PHISPY=1`.

`clean-genome-diffs-gdtools` rewrites `genome-diffs/<sample>.gd` into `genome-diffs-clean/<sample>.gd`,
keeping only mutation entries and writing a dot (`.`) for the id and parent-ids columns, so that editing a
GenomeDiff file under version control does not renumber every entry after the change. It needs a breseq
build that has `gdtools REMOVE --clear-ids`: pass `--breseq-prerelease`, or `--breseq-path DIR` to use a
local build.

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
- `genome_size` — required by autocycler-assemble (can be set per sample in data.csv)
- `BRESEQ_THREADS`, `BRESEQ_OPTIONS`, `NO_DEFAULT_BRESEQ_OPTIONS` — breseq tuning
- `SRA_IGNORE_LAYOUT_MISMATCH` — downgrade the paired/single-end check on `sra://` downloads from an error to a warning

Any of these can also be set per sample in data.csv (see above). Use `--config-add KEY=VALUE` (same as `--config KEY_ADD=VALUE`) to add to data.csv settings instead of overriding them.
