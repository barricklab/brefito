try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "download-data-lftp.smk"

# rule merge_nanopore_reads:
#     input:
#         lambda wildcards: sample_info.get_nanopore_read_list(wildcards.sample) 
#     output:
#         "nanopore_reads_merged/{sample}.fastq.gz"
#     run:
#         # This checks for whether there are two or more, in which case we have to cat
#         if len(input)>1:
#             shell("cat {input} > {output}")
#         else:
#             shell("ln -s ../{input} {output}")

rule trim_nanopore_reads:
    input:
        "nanopore_reads/{sample}.fastq.gz"
    output:
        "nanopore_reads_trimmed/{sample}.fastq.gz"
    log:
        "logs/porechop_{sample}.log"
    conda:
        "../envs/porechop.yml"
    threads: 12
    shell:
        "porechop --discard_middle -i {input} -o {output} > {log} 2>&1"