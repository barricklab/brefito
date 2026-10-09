try: sample_info
except NameError:
    include: "load-sample-info.smk"

# Provides the per-sample merged trimmed reads (nanopore-reads-trimmed-merged/),
# and transitively trim-nanopore-reads.smk / download-data.smk.
include: "filter-nanopore-reads.smk"

rule all_evaluate_nanopore_reads:
    input:
        ["nanopore_read_stats/" + s for s in sample_info.get_samples_with_nanopore_reads()]
    default_target: True

rule evaluate_nanopore_reads:
    input:
        "nanopore-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        dir = directory("nanopore_read_stats/{sample}"),
        file = "nanopore_read_stats/{sample}/NanoStats.txt" 
    log:
        "logs/{sample}/nanoplot.log"
    conda:
        "../envs/nanoplot.yml"
    threads: 1
    shell:
        "NanoPlot -t {threads} --fastq {input} -o {output.dir} --title {wildcards.sample} --tsv_stats  --info_in_report --plots dot --legacy hex --loglength > {log} 2>&1"