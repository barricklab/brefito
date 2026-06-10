try: sample_info
except NameError:
    include: "load-sample-info.smk"

include: "predict-mutations-breseq.smk"

import os.path
from pathlib import Path

# Trees only make sense if all of the genome diffs have been generated relative to the same reference. Abort process if otherwise.
if (sample_info.get_reference_prefix()  != "references"):
    raise ValueError(f"Phylogenetic trees can only be constructed with the same reference!")

PHYLOGENY_REFERENCE=None

if 'PHYLOGENY_REFERENCE' in brefito_config.keys(): # name of the reference file. Assuming it is in `references`
    PHYLOGENY_REFERENCE = brefito_config['PHYLOGENY_REFERENCE']

if PHYLOGENY_REFERENCE is None:
    raise ValueError("Please provide the config value for PHYLOGENY REFERENCE")

if not (Path(f"references/{PHYLOGENY_REFERENCE}")).exists(): # check to see if the file exists
    raise FileNotFoundError(f"File references/{PHYLOGENY_REFERENCE} does not exist.")

PHYLOGENY_OPTIONS=""
if 'PHYLOGENY_OPTIONS' in brefito_config.keys():
    PHYLOGENY_OPTIONS = brefito_config['PHYLOGENY_OPTIONS']

file_prefix = "phylogeny"

default_tree_file_name = "breseq-" + sample_info.get_reference_prefix() + f"/trees/{file_prefix}.tre" # default name for output
default_log_file_name = "logs/gdtools-" + sample_info.get_reference_prefix() + f"{file_prefix}.log" #default name for log files

def get_tree_file_name_prefix(default_tree_file_name): # deciding on the name for the output to avoid overwriting the existing output
    path = Path(default_tree_file_name)

    if not path.exists():
        return file_prefix

    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return f"{file_prefix}_1"
        i += 1

rule generate_phylogenetic_tree:
    input:
        gd = expand("breseq-" + sample_info.get_reference_prefix() + "/gd/{sample}.gd", sample=sample_info.get_sample_list()), # you can't use wildcards here but you can use this expand functionality
        references = f"references/{PHYLOGENY_REFERENCE}"
    output:
        tree_file = "breseq-" + sample_info.get_reference_prefix() + "/trees/" + get_tree_file_name_prefix(default_tree_file_name) + "/" + get_tree_file_name_prefix(default_tree_file_name) + ".tre"
    log:
        "logs/gdtools-" + sample_info.get_reference_prefix() + "-" + get_tree_file_name_prefix(default_tree_file_name) + ".log"
    conda:
        "../envs/breseq-phylip.yml"
    params:
        gd_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/gd"),
        prefix = get_tree_file_name_prefix(default_tree_file_name),
        tree_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/trees/" + get_tree_file_name_prefix(default_tree_file_name) + "/"),
    shell:
        """
        mkdir -p {params.tree_dir}
        gdtools PHYLOGENY -r {input.references} -o {params.prefix} {input.gd}
        mv {params.prefix}.genotypes.txt {params.prefix}.merged.gd {params.prefix}.mutation.key.txt {params.prefix}.phylip.commands.txt {params.prefix}.phylip.output.txt {params.prefix}.sample.key.txt {params.prefix}.tre {params.tree_dir}

        """
