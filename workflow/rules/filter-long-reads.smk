include: "trim-long-reads.smk"

rule gunzip_long_reads:
    input:
        reads = "01_trimmed_reads/{dataset}.long_reads.trimmed.fastq.gz"
    output:
        reads = "01_trimmed_reads/{dataset}.long_reads.trimmed.fastq"
    shell:
        "gunzip {input}"

rule filter_long_reads:
    input:
        "01_trimmed_reads/{sample}.long_reads.trimmed.fastq"
    output:
        "01_trimmed_reads/{sample}.long_reads.trimmed.filtered.fastq"
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input} > {output}"