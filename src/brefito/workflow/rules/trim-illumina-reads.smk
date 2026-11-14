try: sample_info
except NameError: 
    include: "load-sample-info.smk"

try: DOWNLOAD_DATA_INCLUDED
except NameError: 
    include: "download-data.smk"

#supply 'none' to skip
# Picks the rule order, so it can't vary between samples
ILLUMINA_TRIMMING = str(get_uniform_config_value('ILLUMINA_TRIMMING', "fastp"))

# Defaults. FASTP_OPTIONS, FASTP_SE_OPTIONS and FASTP_PE_OPTIONS can be set per sample in
# data.csv. Trimming is per read file, so every sample listing the same read file must
# agree on them.
FASTP_OPTIONS = "--disable_quality_filtering"
FASTP_SE_OPTIONS = ""
FASTP_PE_OPTIONS = "--detect_adapter_for_pe"
if 'FASTP_PE_OPTIONS' in brefito_config.keys():
    print("NOTE: Default fastp PE --detect_adapter_for_pe --disable_quality_filtering options overwritten by user option.")

def get_fastp_option(key, default, read_file_name):
    return get_uniform_config_value(key, default, sample_info.get_samples_for_read_file(read_file_name))



READ_NUMS = ["1", "2"]

if ILLUMINA_TRIMMING.upper()=="NONE":
    print("Illumina reads are not trimmed.")
    ruleorder: trim_illumina_reads_no_trim > trim_PE_illumina_reads_with_fastp
    ruleorder: trim_illumina_reads_no_trim > trim_SE_illumina_reads_with_fastp
elif ILLUMINA_TRIMMING.upper()=="FASTP":
    print("fastp used for trimming illumina reads.")
    # fastp pairs R1 with R2 through a shared {dataset} wildcard. If the names of a pair do
    # not share a base, no fastp job exists for them and snakemake quietly falls back to
    # trim_illumina_reads_no_trim, passing the reads through untrimmed. Check here so that
    # this is an error rather than a silent no-op.
    for this_sample in sample_info.get_sample_list():
        sample_info.get_illumina_PE_read_base_list(this_sample)
    ruleorder: trim_PE_illumina_reads_with_fastp > trim_illumina_reads_no_trim 
    ruleorder: trim_SE_illumina_reads_with_fastp > trim_illumina_reads_no_trim
    print("  Single-end fastp options: " + str(get_config_value('FASTP_SE_OPTIONS', FASTP_SE_OPTIONS)))
    print("  Paired-end fastp options: " + str(get_config_value('FASTP_PE_OPTIONS', FASTP_PE_OPTIONS)))
    print("  Shared     fastp options: " + str(get_config_value('FASTP_OPTIONS', FASTP_OPTIONS)))
    # Check now that samples sharing a read file agree on per-sample fastp options
    for this_sample in sample_info.get_sample_list():
        for this_read_file in sample_info.get_illumina_read_list(this_sample):
            for this_key, this_default in [('FASTP_OPTIONS', FASTP_OPTIONS), ('FASTP_SE_OPTIONS', FASTP_SE_OPTIONS), ('FASTP_PE_OPTIONS', FASTP_PE_OPTIONS)]:
                get_fastp_option(this_key, this_default, this_read_file)
else:
    print("Unknown trimming method requested for --config ILLUMINA_TRIMMING.")
    print("Valid options are: 'none', 'fastp' (default).")
    exit(0)

def all_trimmed_illumina_read_names():
    l = [] 
    for s in sample_info.get_sample_list():
        l = l + ["illumina-reads-trimmed/" + r for r in sample_info.get_illumina_read_list(s) ]
    return(l)

rule trim_all_illimina_reads:
    input:
        all_trimmed_illumina_read_names()

rule trim_PE_illumina_reads_with_fastp:
    input:
        expand("illumina-reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("illumina-reads-trimmed/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/fastp.yml"
    log:
        html = "logs/fastp/{dataset}.RX.html",
        json = "logs/fastp/{dataset}.RX.json",
        log = "logs/fastp/{dataset}.RX.log"
    params:
        fastp_options = lambda wildcards: get_fastp_option('FASTP_OPTIONS', FASTP_OPTIONS, wildcards.dataset + ".R1.fastq.gz"),
        fastp_pe_options = lambda wildcards: get_fastp_option('FASTP_PE_OPTIONS', FASTP_PE_OPTIONS, wildcards.dataset + ".R1.fastq.gz"),
    threads: 12
    shell:
        """
        fastp {params.fastp_options} {params.fastp_pe_options} -j {log.json} -h {log.html} --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} > {log.log} 2>&1
        """

rule trim_SE_illumina_reads_with_fastp:
    input:
        "illumina-reads/{dataset}.SE.fastq.gz"
    output:
        "illumina-reads-trimmed/{dataset}.SE.fastq.gz"
    conda:
        "../envs/fastp.yml"
    log:
        html = "logs/fastp/{dataset}.SE.html",
        json = "logs/fastp/{dataset}.SE.json",
        log = "logs/fastp/{dataset}.SE.log"
    params:
        fastp_options = lambda wildcards: get_fastp_option('FASTP_OPTIONS', FASTP_OPTIONS, wildcards.dataset + ".SE.fastq.gz"),
        fastp_se_options = lambda wildcards: get_fastp_option('FASTP_SE_OPTIONS', FASTP_SE_OPTIONS, wildcards.dataset + ".SE.fastq.gz"),
    threads: 12
    shell:
        """
        fastp {params.fastp_options} {params.fastp_se_options} -j {log.json} -h {log.html} --thread {threads} -i {input} -o {output} > {log.log} 2>&1
        """

rule trim_illumina_reads_no_trim:
    input:
        "illumina-reads/{sample}.fastq.gz"
    output:
        "illumina-reads-trimmed/{sample}.fastq.gz"
    threads: 1
    shell:
        "ln -s ../{input} {output}"