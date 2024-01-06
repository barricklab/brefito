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

rule all:
    input:
        ["breseq/" + s + "/output/output.done" for s in sample_info.get_sample_list()]
    default_target: True

rule predict_mutations_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        references = lambda wildcards: [ "references/" + d for d in sample_info.get_reference_list(wildcards.sample)]
    output:
        breseq = directory("breseq/{sample}"),
        output = directory("output/{sample}"),
        done_file = "breseq/{sample}/output/output.done"
    log: 
        "logs/evaluate_breseq_{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r references/')
    threads: 8
    shell:
        """
        mkdir -p output
        mkdir -p breseq
        breseq -j {threads} {params.automatic_breseq_args} {BRESEQ_OPTIONS} {params.reference_arguments} -o {output.breseq}  {input.reads} > {log} 2>&1
        cp -r {output.breseq}/output {output.output}
        """