include: "align-for-igv.smk"

BREFITO_PACKAGE_PATH = config["brefito_package_path"]

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule evaluate_soft_clipping:
    input:
        bam = "evaluate/aligned_reads/{sample}/{technology}.bam",
        fasta = "evaluate/aligned_reads/{sample}/reference.fasta",
    output:
        csv = "evaluate/soft_clipping/{technology}/{sample}.csv" 
    log:
        "logs/{sample}/evaluate_soft_clipping_{technology}.log"
#    conda:
#        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq SOFT-CLIPPING -o {output.csv} -b {input.bam} -f {input.fasta} > {log} 2>&1
        """

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule summarize_soft_clipping:
    input:
        csv = "evaluate/soft_clipping/{technology}/{sample}.csv" 
    output:
        plot = "evaluate/soft_clipping_summary/{technology}/{sample}_soft_clipping_plot.pdf",
        summary = "evaluate/soft_clipping_summary/{technology}/{sample}_soft_clipping_summary.csv" 
    log:
        "logs/{sample}/summarize_soft_clipping_{technology}.log"
    conda:
        "../envs/rstats.yml"
    threads: 1
    shell:
        """
        Rscript {BREFITO_PACKAGE_PATH}/summarize_soft_clipping.R -i {input.csv} -p {output.plot} -o {output.summary} > {log} 2>&1
        """

