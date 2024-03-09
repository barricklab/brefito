try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-illumina-reads.smk"

import os.path

BRESEQ_OPTIONS = ""
if 'BRESEQ_OPTIONS' in brefito_config.keys():
    BRESEQ_OPTIONS = brefito_config['BRESEQ_OPTIONS'] 

NO_DEFAULT_BRESEQ_OPTIONS = False
if 'NO_DEFAULT_BRESEQ_OPTIONS' in brefito_config.keys():
    NO_DEFAULT_BRESEQ_OPTIONS = bool(brefito_config['NO_DEFAULT_BRESEQ_OPTIONS'])

AUTOMATIC_BRESEQ_ARGS = {}

def find_available_read_files(wildcards):

    illumina_files = [ "illumina-reads-trimmed/" + d for d in sample_info.get_illumina_read_list(wildcards.sample)]

    AUTOMATIC_BRESEQ_ARGS[wildcards.sample] = ""

    return illumina_files

def get_breseq_args(sample):
    return AUTOMATIC_BRESEQ_ARGS[sample]

rule all_predict_mutations_breseq:
    input:
        ["assemblies/" + s + ".fasta.polished" for s in sample_info.get_sample_list()]
    default_target: True

# Add automating apply step
#rule apply_predicted_mutations:


rule apply_mutations_gdtools:
    input:
        gd_file = "breseq-polish" + "/gd/{sample}.gd",
        reference = "assemblies/{sample}.fasta",
    output:
        "assemblies/{sample}.fasta.polished",
    log: 
        "logs/breseq-polish-apply" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        gdtools APPLY -f fasta -r {input.reference} -o  {output} {input.gd_file}
        """

rule predict_mutations_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        reference = "assemblies/{sample}.fasta",
    output:
        breseq_dir = temp(directory("breseq-polish" + "/data/{sample}")),
        html_dir = directory("breseq-polish" + "/html/{sample}"),
        done_file = "breseq-polish" + "/html/{sample}/output.done",
        gd_file = "breseq-polish" + "/gd/{sample}.gd",
        bam = "breseq-polish" + "/data/{sample}/data/reference.bam",
        fasta = "breseq-polish" + "/data/{sample}/data/reference.fasta",
        summary_json = "breseq-polish" + "/data/{sample}/data/summary.json",
    log: 
        "logs/breseq-polish-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        gd_dir = directory("breseq-polish" + "/gd"),
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
    threads: 8
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} {params.automatic_breseq_args} {BRESEQ_OPTIONS} -c {input.reference} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """