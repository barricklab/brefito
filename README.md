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

It's also highly recommended that you set up a common location for snakemake to store conda environments. Otherwise it will re-install each environment in each working directory you use, which is slow and wastes disk space! To do so, set this environment variable in your startup script (`~/.bashrc`, `~/.zshrc`, etc.) to something like this:
```
export SNAKEMAKE_CONDA_PREFIX=$HOME/snakemake_conda_envs
```

## Overview

_brefito_ is a wrapper to make running several related snakemake pipelines easier!

## Working Directory Structure

To get started, you'll need to place your input files in some standard locations withing the main working directory that you will use for running _brefito_.

Input directories of sequencing reads:
* `input/illumina_reads`
* `input/nanopore_reads`

Special input/output directory, depending on the command
* `assemblies`

Input directories for commands that compare genomes:
* `input/references`
* `input/samples`

Output files usually appear in directories of the form `output/*`

## Commands

### `brefito download-reads-lftp`
```
Inputs: data.csv
Outputs: nanopore_reads, illumina_reads
Options: --config bookmark=<lftp bookmark name>
Options: --config bookmark=<lftp bookmark name>
```

### `brefito evaluate-breseq`
```
Inputs:
Output:
Options: --config breseq_options="<breseq_options>"
```