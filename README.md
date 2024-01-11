# brefito

Helper for running snakemake pipelines for assembling, polishing, and annotating bacterial reference sequences and for predicting mutations in a bacterial genome from re-sequencing data.

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

## Specifying Input Data

### Using a CSV file

Each line specifies a sample and a reference or read file or program option associated with that sample.

```
sample,type,setting
Ara-1_75000gen_A,reference,https://raw.githubusercontent.com/barricklab/LTEE/7da91974eafac0c5a8f903ae57275795d4395af2/reference/REL606.gbk
Ara-1_75000gen_A,illumina-R1,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR259/007/SRR2591047/SRR2591047_1.fastq.gz
Ara-1_75000gen_A,illumina-R2,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR259/007/SRR2591047/SRR2591047_2.fastq.gz
```

File locations prefixed with `ftp://` or `http://` will be downloaded using wget. You must have public access to these.

If you have files stored in a private server, you can access them using lftp bookmarks. In this case, us the prefix `lftp@[bookmark_name]://`.

Currently, you can use the `*` wildcard name for sample when you want a file to be associated with all samples, BUT this must be at the end of the file. (It only gets applied to samples it already knows about when it reads that line.)

See the Examples directory for `data.csv` files you can use for testing and as templates.

### Using a Directory Structure (CURRENTLY BROKEN)

Instead of creating the CSV file, you can place input files in some standard locations within the main working directory that you will use for running _brefito_. This option is more limited in terms of flexibility. For example, you can't use multiple reference files or read file sets for a sample, you'll need to make just one per sample.

Input directories of sequencing reads:
* `input/illumina_reads`
* `input/nanopore_reads`

Special input/output directory, depending on the command
* `assemblies`

Input directories for commands that compare genomes:
* `input/references`
* `input/samples`

Global reference file for comparisons
* `reference.fasta`

Output files usually appear in directories of the form `output/*`

## Commands

### `brefito download-reads-lftp`

Downloads the files specified in a sample

```
Inputs: data.csv
Outputs: nanopore_reads, illumina_reads
Options: --config data_csv=<path> (Default is data.csv)
         --resources connections=<int>
```

The `connections` resource controls how many simultaneous download jobs will be used. By default it is 1. Be careful to not make it too high and overload your system.

### `brefito predict-mutations-breseq`

Runs `breseq` using the reference files and trimmed read files.

```
Inputs: data.csv
Outputs: breseq_output, output
Options: <same options as download-reads-lftp>
         --config breseq_options="<breseq_options>"
         --config no_default_breseq_options=<bool>
        

```