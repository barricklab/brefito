try: sample_info
except NameError:
    include: "load-sample-info.smk"

# Provides the per-sample merged trimmed reads (nanopore-reads-trimmed-merged/),
# and transitively trim-nanopore-reads.smk / download-data.smk.
include: "filter-nanopore-reads.smk"

# NanoPlot summarizes any long-read FASTQ, so PacBio reads are evaluated the same way, into
# <type>_read_stats/ (nanopore_read_stats/ for nanopore reads, as before).
rule all_evaluate_nanopore_reads:
    input:
        [sample_info.get_long_read_type(s) + "_read_stats/" + s for s in sample_info.get_samples_with_long_reads()]
    default_target: True

rule evaluate_nanopore_reads:
    input:
        "{lrt}-reads-trimmed-merged/{sample}.fastq.gz"
    output:
        dir = directory("{lrt}_read_stats/{sample}"),
        file = "{lrt}_read_stats/{sample}/NanoStats.txt"
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    log:
        "logs/{sample}/nanoplot-{lrt}.log"
    conda:
        "../envs/nanoplot.yml"
    threads: 1
    shell:
        "NanoPlot -t {threads} --fastq {input} -o {output.dir} --title {wildcards.sample} --tsv_stats  --info_in_report --plots dot --legacy hex --loglength > {log} 2>&1"