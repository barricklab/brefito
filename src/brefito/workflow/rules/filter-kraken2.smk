try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

# Different options

if 'KRAKEN2_FILTER_TAXID' in brefito_config.keys():
    KRAKEN2_FILTER_TAXID = brefito_config['KRAKEN2_FILTER_TAXID']
else:
    raise Exception("Must provide --config KRAKEN2_FILTER_TAXID=<NCBI tax id>")

if 'KRAKEN2_DB' in brefito_config.keys():
    KRAKEN2_DB = brefito_config['KRAKEN2_DB']
else:
    raise Exception("Must provide --config KRAKEN2_DB=<path to database>")

rule all_filter_taxid_kraken2:
    input:
        lambda wildcards: ["filtered-fastq/" + s + ".fastq.gz" for s in sample_info.sample_list]
    default_target: KRAKEN2_FILTER_TAXID == ""
    conda:
        "../envs/kraken2.yml"
    shell:
        "k2 clean --stop-daemon"

rule gzip_filtered_fastq:
    input:
        "filter-kraken2/filtered-fastq/{sample}.fastq"
    output:
        "filtered-fastq/{sample}.fastq.gz"
    conda:
        "../envs/download.yml"
    threads: 8
    shell:
        "pigz -p {threads} -c {input} > {output}"


rule filter_taxid_krakentools:
    input:
        fastq = lambda wildcards: ["illumina-reads-trimmed/" + r for r in sample_info.get_illumina_read_list(wildcards.sample)],
        kraken_output = "filter-kraken2/kracken2_output/{sample}.txt",
        kraken_report = "filter-kraken2/kracken2_report/{sample}.txt"
    output:
        temp("filter-kraken2/filtered-fastq/{sample}.fastq")
    log:
        "logs/filter-taxid-krakentools-{sample}.log"
    params:
        fastq_arg = lambda wildcards,input: [ "-s " + r for r in input.fastq]
    conda:
        "../envs/kraken2.yml"
    shell:
       "extract_kraken_reads.py {params.fastq_arg} -o {output} -r {input.kraken_report} -k {input.kraken_output} -t {KRAKEN2_FILTER_TAXID} --include-children --fastq-output"

rule filter_kraken2:
    input:
        lambda wildcards: ["illumina-reads-trimmed/" + r for r in sample_info.get_illumina_read_list(wildcards.sample)]
    output:
        report = "filter-kraken2/kracken2_report/{sample}.txt",
        output = "filter-kraken2/kracken2_output/{sample}.txt"
    log:
        "logs/filter-kraken2-{sample}.log"
    conda:
        "../envs/kraken2.yml"
    threads: 6
    resources:
        singleton=1
    shell:
        "k2 classify --use-daemon --db {KRAKEN2_DB} --threads {threads} --report {output.report} --output {output.output}  {input} > {log} 2>&1"
