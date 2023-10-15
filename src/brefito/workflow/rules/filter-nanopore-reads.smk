include: "trim-nanopore-reads.smk"

#rule gunzip_trimmed_nanopore_reads:
#    input:
#        "01_trimmed_nanopore_reads/{dataset}.fastq.gz"
#    output:
#        temp("01_trimmed_nanopore_reads/{dataset}.fastq")
#    shell:
#        "gunzip < {input} > {output}"


rule filter_nanopore_reads:
    input:
        "01_trimmed_nanopore_reads/{sample}.fastq.gz"
    output:
        "02_filtered_nanopore_reads/{sample}.fastq"
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input} > {output}"