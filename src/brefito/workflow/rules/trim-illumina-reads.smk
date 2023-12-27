READ_NUMS = ["1", "2"]

rule trim_PE_illumina_reads_with_fastp:
    input:
        expand("illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("illumina_reads_trimmed/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/fastp.yml"
    log:
        html = "logs/fastp/{dataset}.RX.html",
        json = "logs/fastp/{dataset}.RX.json",
        log = "logs/fastp/{dataset}.RX.log"
    threads: 12
    shell:
        """
        fastp --detect_adapter_for_pe -j {log.json} -h {log.html} --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} > {log.log} 2>&1
        """

rule trim_SE_illumina_reads_with_fastp:
    input:
        "illumina_reads/{dataset}.SE.fastq.gz"
    output:
        "illumina_reads_trimmed/{dataset}.SE.fastq.gz"
    conda:
        "../envs/fastp.yml"
    log:
        html = "logs/fastp/{dataset}.SE.html",
        json = "logs/fastp/{dataset}.SE.json",
        log = "logs/fastp/{dataset}.SE.log"
    threads: 12
    shell:
        """
        fastp -j {log.json} -h {log.html} --thread {threads} -i {input} -o {output} > {log.log} 2>&1
        """
