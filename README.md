# brefito

Helper for running snakemake pipelines for assembling, polishing, checking, and annotating bacterial reference sequences and for predicting mutations in a bacterial genome from re-sequencing data.

## Installation

Clone the repository and change to the main directory
```
git clone https://github.com/barricklab/brefito.git
cd brefito
```

Install conda and mamba.

Then, create an environment for _brefito_ thet includes its requirements and install it in this environment.

```
mamba env create -f environment.yml
conda activate brefito
pip install .
```

Now you should be able to run the _brefito_ script and see the help like this!
```
brefito
```

It's also *highly recommended* that you set up a common location for snakemake to store conda environments. Otherwise it will re-install each environment in each working directory you use, which is slow and wastes disk space! To do so, set this environment variable in your startup script (`~/.bashrc`, `~/.zshrc`, etc.) to something like this:
```
export SNAKEMAKE_CONDA_PREFIX=$HOME/snakemake_conda_envs
```

## Overview

_brefito_ is a wrapper to make running several related snakemake pipelines easier!

## Usage

```
usage: brefito command
```

_brefito_ has several options that are pass-throughs to snakemake.

```
  --config CONFIG
  --resources RESOURCES
  --rerun-incomplete
  --unlock
  --keep-going
  --dry-run
```
These can allow you to change certain resources or configuration variables (described below for each workflow).

They can also allow you to resume and control execution after failed/interrupted runs.

## Specifying Input Data

### Using a CSV file

Each line specifies a sample and a reference or read file or program option associated with that sample.

```
sample,type,setting
sample,type,setting
Ara-1_50000gen_11331,reference,https://raw.githubusercontent.com/barricklab/LTEE/7da91974eafac0c5a8f903ae57275795d4395af2/reference/REL606.gbk
Ara-1_50000gen_11331,illumina-R1,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR259/007/SRR2591047/SRR2591047_1.fastq.gz
Ara-1_50000gen_11331,illumina-R2,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR259/007/SRR2591047/SRR2591047_2.fastq.gz
```

File locations prefixed with `ftp://` or `http://` will be downloaded using wget. You must have public access to these.

If you have files stored in a private server, you can access them using lftp bookmarks. In this case, us the prefix `lftp@[bookmark_name]://`.

Currently, you can use the `*` wildcard name for sample when you want a file to be associated with all samples, BUT this must be at the end of the file. (It only gets applied to samples it already knows about when it reads that line.)

See the Examples directory for `data.csv` files you can use for testing and as templates.

## Commands

### `brefito download-data`

Downloads the files specified in a sample

```
Inputs: data.csv
Outputs: nanopore_reads, illumina_reads
Options: --config data_csv=<path> (Default is data.csv)
         --resources connections=<int>
```

The `connections` resource controls how many simultaneous download jobs will be used. By default it is 1. Be careful to not make it too high and overload your system!

### `brefito merge-reads`

Merges the raw reads corresponding to each sample into one file per type of read.

```
Inputs: data.csv
Outputs: merged_reads
Options: <same options as download-data>
         --config breseq_options="<breseq_options>"
```

### `brefito merge-trimmed-reads`

Merges the trimmed reads corresponding to each sample into one file per type of read.

```
Inputs: data.csv
Outputs: merged_reads_trimmed
Options: <same options as download-data>
         --config breseq_options="<breseq_options>"
```

### `brefito predict-mutations-breseq`

Runs `breseq` using the reference files and trimmed read files.

```
Inputs: data.csv
Outputs: breseq_reference/data, breseq_reference/html, breseq_reference/gd 
Options: <same options as download-data>
         --config breseq_options="<breseq_options>
            Options that get passed to breseq
         --config no_default_breseq_options=<bool>
            Don't pass the default option of -x to breseq when using nanopore reads
```

### `brefito coverage-plots-breseq`

Runs `breseq BAM2COV` to create coverage plots tiling the reference genome.

```
inputs: breseq_reference/data, breseq_reference/html, breseq_reference/gd
Outputs: breseq_reference/cov
Options: <same options as download-data>
         --config breseq_options="<breseq_options>
            Options that get passed to breseq
         --config no_default_breseq_options=<bool>
            Don't pass the default option of -x to breseq when using nanopore reads
```

### `brefito align-reads`

Generates files that can be loaded in IGV to view sequences (FASTA/FAI), reads (BAM/BAI) and annotations (GFF). Runs minimap2 for nanopore reads and bowtie2 for illumina reads for mapping to the provided reference.

```
Inputs: data.csv
Outputs: align_reads_reference/data
Options: <same options as download-data>
```

### `brefito check-soft-clipping`

Analyzes and plots soft-clipped reads after mapping.

```
Inputs: align_reads_reference/data
Outputs: align_reads_reference/soft-clipping
Options: <same options as download-data>
```

### `brefito mutate-genomes-gdtools`

Uses `gdtools` from _breseq_ to apply the GenomeDiff files in `genome_diff` to generate updated reference genomes that include those mutations. One GenomeDiff file is expected per sample with the `*.gd` file ending. These could be copied from a `breseq-*/gd` directory and then manually edited to curate the mutations they describe.

```
Inputs: data.csv, genome_diffs/*.gd
Outputs: mutants
Options: <same options as download-data>
```

### `brefito annotate-genomes`

Combines annotations of genes from `prokka` with annotations of IS elements from `isescan` into a final Genbank file for each sample.

```
Inputs: data.csv
Outputs: annotated_references
Options: <same options as download-data>
```

### Using a different directory of reference sequences

Many workflows that operate, by default, on the reference sequences in the `references` directory can be used on a different directory of reference files. The main cases when this is useful is when re-running analyses on `assemblies` or `mutants` generated by other workflows to check them for accuracy.

In this case you can use the normal command with another dash and the desired directory of reference files appended. For example: `predict-mutations-breseq-assemblies` or `predict-mutations-breseq-mutants`, `coverage-plots-breseq-assemblies` or `coverage-plots-breseq-mutants`, `align-reads-assemblies`, `check-soft-clipping-assemblies`, etc.

These commands will also work on input directories other than `assemblies` or `mutants`, just replace the same part of the command with your input folder name.

For example, you can "chain" execution of `annotate-genomes-assemblies` and `predict-mutations-breseq-annotated_assemblies` to use annotated assemblies as the reference sequences for running _breseq_.
