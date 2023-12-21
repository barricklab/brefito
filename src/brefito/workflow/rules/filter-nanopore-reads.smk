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

rule filter_nanopore_reads:
    input:
        "01_trimmed_nanopore_reads/{sample}.fastq.gz"
    output:
        "02_filtered_nanopore_reads/{sample}.fastq"
    log:
        "logs/filtlong_{sample}.fastq"
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent {FILTLONG_KEEP_PERCENT} --target_bases {FILTLONG_TARGET_BASES} {input} > {output} 2> {log}"