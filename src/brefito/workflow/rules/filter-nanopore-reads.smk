try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"

do_filtering = True
# Picks the rule order, so it can't vary between samples
NANOPORE_FILTERING = str(get_uniform_config_value('NANOPORE_FILTERING', "TRUE"))

if NANOPORE_FILTERING.upper() == "NONE" or NANOPORE_FILTERING.upper() == "FALSE" or NANOPORE_FILTERING.upper() == "F" or NANOPORE_FILTERING == "0":
    print("Nanopore reads are not filtered.")
    do_filtering = False
    ruleorder: filter_nanopore_reads_no_filter > filter_nanopore_reads
else:
    ruleorder: filter_nanopore_reads > filter_nanopore_reads_no_filter

FILTLONG_OPTIONS = "--min_length 1000 --keep_percent 90"
if (do_filtering):
    if "FILTLONG_OPTIONS" not in brefito_config.keys():
        print("Filtering nanopore reads using default filtlong command-line options (unless set per sample in data.csv):")
        print("  " + FILTLONG_OPTIONS)
        print("To change these options add something like this to your brefito command")
        print("  --config FILTLONG_OPTIONS=\"--min_length 5000 --keep_percent 95 --target_bases 500000000\"")
    else:
        print("Filtering nanopore reads with user-specified filtlong command-line options:")
        print("  " + brefito_config["FILTLONG_OPTIONS"])

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
    params:
        filtlong_options = lambda wildcards: get_config_value('FILTLONG_OPTIONS', FILTLONG_OPTIONS, wildcards.sample)
    shell:
        "filtlong {params.filtlong_options} {input} > {output} 2> {log}"


rule filter_nanopore_reads_no_filter:
    input:
        # the same trimmed and merged reads filtlong would have filtered
        "nanopore-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        "nanopore-reads-filtered/{sample}.fastq"
    threads: 1
    shell:
        "gunzip -c {input} > {output}"