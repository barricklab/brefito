READ_NUMS = ["1", "2"]

rule trim_with_fastp:
    input:
        expand("illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("01_trimmed_illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/fastp.yml"
    log:
        html = "logs/fastp/{dataset}.html",
        json = "logs/fastp/{dataset}.json",
        log = "logs/fastp/{dataset}.log"
    threads: 12
    shell:
        """
        fastp --detect_adapter_for_pe -j {log.json} -h {log.html} --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} > {log.log} 2>&1
        """
