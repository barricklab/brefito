include: "predict-mutations-breseq.smk"

rule all_coverage_plots_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/cov/" + s + "/coverage.done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to create  output files for each reference by default if no argument provided...
rule coverage_plots_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        breseq_summary_json = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json"
    output:
        dir = directory("breseq-" + sample_info.get_reference_prefix() + "/cov/{sample}"),
        #whole_genome_file = "evaluate/coverage_plots/{technology}/{sample}/all.pdf",
        done_file = "breseq-" + sample_info.get_reference_prefix() + "/cov/{sample}/coverage.done"
    log:
        "logs/coverage-plots-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq bam2cov -o {output.dir} -x tiling_ -b {input.bam} -f {input.fasta} -a --tile-size 10000 --tile-overlap 2000 -p 1000 -j {input.breseq_summary_json} > {log} 2>&1
        breseq bam2cov -o {output.dir} -x summary_ -b {input.bam} -f {input.fasta} -a -j {input.breseq_summary_json} -p 10000 >> {log} 2>&1
        touch {output.done_file}
        """
