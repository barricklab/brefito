try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

import os.path

BRESEQ_THREADS = int(get_uniform_config_value('BRESEQ_THREADS', 8))

# BRESEQ_OPTIONS can be set per sample in data.csv
def get_breseq_options(sample):
    return get_config_value('BRESEQ_OPTIONS', "", sample)

def find_available_read_files(wildcards):

    # long reads of any type; breseq aligns them itself here (for PacBio CLR reads, the
    # minimap2 route in predict-mutations-minimap2-breseq usually maps them better)
    long_read_type = sample_info.get_long_read_type(wildcards.sample)
    long_read_files = [ long_read_type + "-reads-trimmed/" + d for d in sample_info.get_long_read_list(wildcards.sample)] if long_read_type else []
    illumina_files = [ "illumina-reads-trimmed/" + d for d in sample_info.get_illumina_read_list(wildcards.sample)]

    return illumina_files + long_read_files

def get_breseq_args(sample):
    return ""

rule all_predict_mutations_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/html/" + s + "/output.done" for s in sample_info.get_sample_list()]
    default_target: True

rule predict_mutations_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        breseq_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/data/{sample}"),
        done_file = "breseq-" + sample_info.get_reference_prefix() + "/html/{sample}/output.done",
        gd_file = "breseq-" + sample_info.get_reference_prefix() + "/gd/{sample}.gd",
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        summary_json = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json",
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        html_dir = lambda wildcards: directory("breseq-" + sample_info.get_reference_prefix() + "/html/" + wildcards.sample),
        gd_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/gd"),
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
        breseq_options = lambda wildcards: get_breseq_options(wildcards.sample),
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: BRESEQ_THREADS
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} {params.automatic_breseq_args} {params.breseq_options} {params.reference_arguments} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {params.html_dir}
        mv {output.breseq_dir}/output {params.html_dir}
        """