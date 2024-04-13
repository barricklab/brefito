try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"

#This is 100x nominal coverage of a 5 Mbp genome
FILTLONG_TARGET_BASES = 500000000
if "filtlong_target_bases" in config.keys():
    FILTLONG_TARGET_BASES = config["filtlong_target_bases"]

FILTLONG_KEEP_PERCENT = 90
if "filtlong_keep_percent" in config.keys():
    FILTLONG_KEEP_PERCENT = config["filtlong_keep_percent"]

rule all_filter_nanopore_reads:
    input:
        ["nanopore-reads-filtered/" + s + ".fastq" for s in sample_info.get_sample_list()]
    default_target: True

rule merge_nanopore_reads_for_filtering:
    input:
        lambda wildcards: ["nanopore-reads-trimmed/" + r  for r in sample_info.get_nanopore_read_list(wildcards.sample)]
    output:
        temp("nanopore-reads-trimmed-merged/{sample}.fastq.gz")
    log:
        "logs/merge_nanopore_reads_for_filtering_{sample}.log"
    shell:
        "cat {input} > {output} 2> {log}"

rule filter_nanopore_reads:
    input:
        "nanopore-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        "nanopore-reads-filtered/{sample}.fastq"
    log:
        "logs/filtlong_{sample}.log"
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent {FILTLONG_KEEP_PERCENT} --target_bases {FILTLONG_TARGET_BASES} {input} > {output} 2> {log}"