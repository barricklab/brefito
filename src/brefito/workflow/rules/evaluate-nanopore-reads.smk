rule evaluate_nanopore_reads:
    input:
        "01_trimmed_nanopore_reads/{sample}.fastq"
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