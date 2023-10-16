include: "align-for-igv.smk"
include: "evaluate-breseq.smk"

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule evaluate_coverage:
    input:
        bam = "evaluate/aligned_reads/{sample}/{technology}.bam",
        fasta = "evaluate/aligned_reads/{sample}/reference.fasta",
        breseq_summary_json = "evaluate/breseq_output/{sample}_{technology}_reads/summary.json"
    output:
        dir = directory("evaluate/coverage_plots/{technology}/{sample}"),
        done_file = "evaluate/coverage_plots/{technology}/{sample}/coverage.done" 
    log:
        "logs/{sample}/evaluate_coverage_{technology}.log"
#    Official breseq doesn't support plotting entire sequence (yet)
    conda:
        "../envs/rstats.yml"
    threads: 1
    shell:
        """
        mkdir -p {output.dir}
        # tiles across sequences
        breseq bam2cov -o {output.dir} -b {input.bam} -f {input.fasta} -a --tile-size 10000 --tile-overlap 2000 -p 1000 -j {input.breseq_summary_json} > {log} 2>&1
        # entire sequences
        breseq bam2cov -o {output.dir} -b {input.bam} -f {input.fasta} -a -p 10000 -j {input.breseq_summary_json} > {log} 2>&1
        touch {output.done_file}
        """
        #breseq bam2cov -o {output.whole_genome_file} -b {input.bam} -f {input.fasta} -a -j {input.breseq_summary_json} -p 10000 > {log} 2>&1
