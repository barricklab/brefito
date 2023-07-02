READ_NUMS = ["1", "2"]

rule trim_with_fastp:
    input:
        expand("illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("01_trimmed_illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/fastp.yml"
    threads: 12
    shell:
        "fastp --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]}"
