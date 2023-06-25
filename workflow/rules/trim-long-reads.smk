rule trim_long_reads:
    input:
        "input/{sample}.long_reads.fastq.gz"
    output:
        "01_trimmed_reads/{sample}.long_reads.trimmed.fastq.gz"
    conda:
        "../envs/porechop.yml"
    threads: 16
    shell:
        "porechop --discard_middle -i {input} -o {output}"