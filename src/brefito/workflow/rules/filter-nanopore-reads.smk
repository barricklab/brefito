try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"

FILTLONG_OPTIONS = "--min_length 1000 --keep_percent 90"
if "FILTLONG_OPTIONS" not in config.keys():
    print("Filtering nanopore reads using default filtlong command-line options:")
    print("  " + FILTLONG_OPTIONS)
    print("To change these options add something like this to your brefito command")
    print("  --config FILTLONG_OPTIONS=\"--min_length 5000 --keep_percent 95 --target_bases 500000000\"")
else:
    FILTLONG_OPTIONS = config["FILTLONG_OPTIONS"]
    print("Filtering nanopore reads with user-specified filtlong command-line options:")
    print("  " + FILTLONG_OPTIONS)


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
        "filtlong {FILTLONG_OPTIONS} {input} > {output} 2> {log}"