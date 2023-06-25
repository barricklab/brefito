include: "trim-long-reads.smk"

rule filter_long_reads:
    input:
        "01_trimmed_reads/{sample}.long_reads.trimmed.fastq.gz"
    output:
        "01_trimmed_reads/{sample}.long_reads.trimmed.filtered.fastq"
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input} > {output}"