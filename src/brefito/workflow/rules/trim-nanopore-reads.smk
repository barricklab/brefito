try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "download-data.smk"

def get_all_trimmed_nanopore_read_names():
    nanopore_files = []
    for sample in sample_info.get_sample_list():
        nanopore_files = nanopore_files + [ "nanopore-reads-trimmed/" + d for d in sample_info.get_nanopore_read_list(sample)] 
    return nanopore_files

rule all_trim_nanopore_reads:
    input:
        get_all_trimmed_nanopore_read_names()
    default_target: True

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