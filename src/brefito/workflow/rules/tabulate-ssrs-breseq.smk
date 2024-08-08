try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "predict-mutations-breseq.smk"

SSR_MINIMUM_LENGTH = 6
if 'SSR_MINIMUM_LENGTH' in brefito_config.keys():
    SSR_MINIMUM_LENGTH = brefito_config['SSR_MINIMUM_LENGTH'] 

import os.path

rule all_tabulate_ssrs_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/ssrs/" + s + ".csv" for s in sample_info.get_sample_list()]
    default_target: True

rule tabulate_ssrs_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "breseq-" + sample_info.get_reference_prefix() + "/ssrs/{sample}.csv"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}-tabulate-ssrs.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq CL-TABULATE --strict -m {SSR_MINIMUM_LENGTH} -o {output} --bam {input.bam} --fasta {input.fasta} {params.reference_arguments} > {log} 2>&1
        """
