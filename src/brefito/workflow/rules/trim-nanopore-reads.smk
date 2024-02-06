try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "download-data.smk"

rule trim_nanopore_reads:
    input:
        "nanopore-reads/{sample}.fastq.gz"
    output:
        "nanopore-reads-trimmed/{sample}.fastq.gz"
    log:
        "logs/porechop-{sample}.log"
    conda:
        "../envs/porechop.yml"
    threads: 12
    shell:
        "porechop --discard_middle -i {input} -o {output} > {log} 2>&1"