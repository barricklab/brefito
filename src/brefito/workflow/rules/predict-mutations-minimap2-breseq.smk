try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "align-reads.smk"

import os.path

BRESEQ_THREADS = int(get_uniform_config_value('BRESEQ_THREADS', 8))

# BRESEQ_OPTIONS and NO_DEFAULT_BRESEQ_OPTIONS can be set per sample in data.csv
def get_breseq_options(sample):
    return get_config_value('BRESEQ_OPTIONS', "", sample)

def find_available_read_files(wildcards):

    temp("aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/nanopore_reads.{reads}.sam")

    nanopore_files = [ "aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + wildcards.sample + "/nanopore_reads." + d + ".sam" for d in sample_info.get_nanopore_read_base_list(wildcards.sample)] 

    return nanopore_files

def get_breseq_args(sample):
    if config_is_true("NO_DEFAULT_BRESEQ_OPTIONS", sample=sample):
        return ""
    return "-x --require-match-fraction 0.0"

rule all_predict_mutations_minimap2_breseq:
    input:
        ["minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/" + s + "/output.done" for s in sample_info.get_sample_list()]
    default_target: True

rule predict_mutations_minimap2_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        breseq_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}"),
        html_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/{sample}"),
        done_file = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/{sample}/output.done",
        gd_file = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/gd/{sample}.gd",
        bam = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        summary_json = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json",
    log: 
        "logs/minimap2-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        gd_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/gd"),
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
        breseq_options = lambda wildcards: get_breseq_options(wildcards.sample),
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: BRESEQ_THREADS
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} --aligned-sam {params.automatic_breseq_args} {params.breseq_options} {params.reference_arguments} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """