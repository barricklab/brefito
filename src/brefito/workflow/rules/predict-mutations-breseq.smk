try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

import os.path

BRESEQ_OPTIONS = ""
if 'breseq_options' in config.keys():
    BRESEQ_OPTIONS = config['breseq_options'] 

NO_DEFAULT_BRESEQ_OPTIONS = False
if 'no_default_breseq_options' in config.keys():
    NO_DEFAULT_BRESEQ_OPTIONS = bool(config['no_default_breseq_options'])

AUTOMATIC_BRESEQ_ARGS = {}

def find_available_read_files(wildcards):

    nanopore_files = [ "nanopore_reads_trimmed/" + d for d in sample_info.get_nanopore_read_list(wildcards.sample)]
    #nanopore_files = [os.path.join("nanopore_reads_trimmed", os.path.basename(d)) for d in sample_info.get_nanopore_read_list(wildcards.sample)]
 
    illumina_files = [ "illumina_reads_trimmed/" + d for d in sample_info.get_illumina_read_list(wildcards.sample)]
    #illumina_files = [os.path.join("illumina_reads_trimmed", os.path.basename(d)) for d in sample_info.get_illumina_read_list(wildcards.sample)]

    AUTOMATIC_BRESEQ_ARGS[wildcards.sample] = "";

    if not NO_DEFAULT_BRESEQ_OPTIONS:
        if (len(nanopore_files) > 0):
            AUTOMATIC_BRESEQ_ARGS[wildcards.sample] = AUTOMATIC_BRESEQ_ARGS[wildcards.sample] + "-x"

    #print(illumina_files)
    #print(nanopore_files)

    return illumina_files + nanopore_files

def get_breseq_args(sample):
    #print(AUTOMATIC_BRESEQ_ARGS[sample])
    return AUTOMATIC_BRESEQ_ARGS[sample]

rule all_predict_mutations_breseq:
    input:
        ["breseq_" + sample_info.get_reference_prefix() + "/html/" + s + "/output.done" for s in sample_info.get_sample_list()]
    default_target: True

rule predict_mutations_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        breseq_dir = directory("breseq_" + sample_info.get_reference_prefix() + "/data/{sample}"),
        html_dir = directory("breseq_" + sample_info.get_reference_prefix() + "/html/{sample}"),
        done_file = "breseq_" + sample_info.get_reference_prefix() + "/html/{sample}/output.done",
        gd_file = "breseq_" + sample_info.get_reference_prefix() + "/gd/{sample}.gd"
    log: 
        "logs/breseq_" + sample_info.get_reference_prefix() + "_{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        gd_dir = directory("breseq_" + sample_info.get_reference_prefix() + "/gd"),
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: 8
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} {params.automatic_breseq_args} {BRESEQ_OPTIONS} {params.reference_arguments} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """