try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"

# Filters long reads of every type (LONG_READ_TYPES), each in its own
# <type>-reads-filtered/ directory. NANOPORE_FILTERING turns filtering on or off for all of them.
do_filtering = True
# Picks the rule order, so it can't vary between samples
NANOPORE_FILTERING = str(get_uniform_config_value('NANOPORE_FILTERING', "TRUE"))

if NANOPORE_FILTERING.upper() == "NONE" or NANOPORE_FILTERING.upper() == "FALSE" or NANOPORE_FILTERING.upper() == "F" or NANOPORE_FILTERING == "0":
    print("Long reads are not filtered.")
    do_filtering = False
    ruleorder: filter_nanopore_reads_no_filter > filter_nanopore_reads
else:
    ruleorder: filter_nanopore_reads > filter_nanopore_reads_no_filter

FILTLONG_OPTIONS = "--min_length 1000 --keep_percent 90"
# HiFi reads are uniformly accurate, so ranking them by quality (--keep_percent) only throws away
# good reads; filter them on length alone.
FILTLONG_OPTIONS_HIFI = "--min_length 1000"

def get_filtlong_options(sample):
    default = FILTLONG_OPTIONS_HIFI if sample_info.get_long_read_type(sample) == "pacbio-hifi" else FILTLONG_OPTIONS
    return get_config_value('FILTLONG_OPTIONS', default, sample)

if (do_filtering):
    if "FILTLONG_OPTIONS" not in brefito_config.keys():
        print("Filtering long reads using default filtlong command-line options (unless set per sample in data.csv):")
        print("  " + FILTLONG_OPTIONS + "   (PacBio HiFi: " + FILTLONG_OPTIONS_HIFI + ")")
        print("To change these options add something like this to your brefito command")
        print("  --config FILTLONG_OPTIONS=\"--min_length 5000 --keep_percent 95 --target_bases 500000000\"")
    else:
        print("Filtering long reads with user-specified filtlong command-line options:")
        print("  " + brefito_config["FILTLONG_OPTIONS"])

# Only samples that have long reads (a sample with Illumina reads only has nothing to filter)
rule all_filter_nanopore_reads:
    input:
        [sample_info.get_long_read_path(s, "filtered") for s in sample_info.get_samples_with_long_reads()]
    default_target: True

rule merge_nanopore_reads_for_filtering:
    input:
        lambda wildcards: [wildcards.lrt + "-reads-trimmed/" + r for r in sample_info.get_long_read_list(wildcards.sample)]
    output:
        temp("{lrt}-reads-trimmed-merged/{sample}.fastq.gz")
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    log:
        "logs/merge_{lrt}_reads_for_filtering_{sample}.log"
    shell:
        "cat {input} > {output} 2> {log}"

rule filter_nanopore_reads:
    input:
        "{lrt}-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        "{lrt}-reads-filtered/{sample}.fastq"
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    log:
        "logs/filtlong_{lrt}_{sample}.log"
    conda:
        "../envs/filtlong.yml"
    params:
        filtlong_options = lambda wildcards: get_filtlong_options(wildcards.sample)
    shell:
        "filtlong {params.filtlong_options} {input} > {output} 2> {log}"


rule filter_nanopore_reads_no_filter:
    input:
        # the same trimmed and merged reads filtlong would have filtered
        "{lrt}-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        "{lrt}-reads-filtered/{sample}.fastq"
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    threads: 1
    shell:
        "gunzip -c {input} > {output}"