# brefito

Snakemake pipeline for assembling, polishing, and annotating bacterial reference sequences.

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
pip install -e .
```

Now you should be able to run the _brefito_ script and see the help like this!
```
brefito
```

## Overview

_brefito_ is a wrapper to make running several related snakemake pipelines easier!