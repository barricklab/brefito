rule trim_long_reads:
    input:
        "nanopore_reads/{sample}.fastq.gz"
    output:
        "01_trimmed_nanopore_reads/{sample}.fastq.gz"
    conda:
        "../envs/porechop.yml"
    threads: 16
    shell:
        "porechop --discard_middle -i {input} -o {output}"