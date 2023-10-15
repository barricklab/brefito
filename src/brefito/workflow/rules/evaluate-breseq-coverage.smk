include: "evaluate-breseq.smk"

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule evaluate_breseq_coverage:
    input:
        bam = "evaluate/breseq_data/{sample}_{technology}_reads/reference.bam",
        fasta = "evaluate/breseq_data/{sample}_{technology}_reads/reference.fasta",
        breseq_summary_json = "evaluate/breseq_output/{sample}_{technology}_reads/summary.json"
    output:
        dir = directory("evaluate/coverage_plots/breseq_{technology}/{sample}"),
        #whole_genome_file = "evaluate/coverage_plots/{technology}/{sample}/all.pdf",
        done_file = "evaluate/coverage_plots/breseq_{technology}/{sample}/coverage.done" 
    log:
        "logs/{sample}/evaluate_breseq_coverage_{technology}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        mkdir -p {output.dir}
        breseq bam2cov -o {output.dir} -b {input.bam} -f {input.fasta} -a --tile-size 10000 --tile-overlap 2000 -p 1000 -j {input.breseq_summary_json} > {log} 2>&1
        touch {output.done_file}
        """
        #breseq bam2cov -o {output.whole_genome_file} -b {input.bam} -f {input.fasta} -a -j {input.breseq_summary_json} -p 10000 > {log} 2>&1
