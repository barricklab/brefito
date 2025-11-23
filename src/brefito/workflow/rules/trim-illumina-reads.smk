try: sample_info
except NameError: 
    include: "load-sample-info.smk"

try: DOWNLOAD_DATA_INCLUDED
except NameError: 
    include: "download-data.smk"

READ_NUMS = ["1", "2"]

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
    threads: 12
    shell:
        """
        fastp --detect_adapter_for_pe -j {log.json} -h {log.html} --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} > {log.log} 2>&1
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
    threads: 12
    shell:
        """
        fastp -j {log.json} -h {log.html} --thread {threads} -i {input} -o {output} > {log.log} 2>&1
        """
