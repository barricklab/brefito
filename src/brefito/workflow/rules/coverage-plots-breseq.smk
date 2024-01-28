include: "predict-mutations-breseq.smk"

rule all_coverage_plots_breseq:
    input:
        ["breseq_" + sample_info.get_reference_prefix() + "/cov/" + s + "/coverage.done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule coverage_plots_breseq:
    input:
        bam = "breseq_" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq_" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        breseq_summary_json = "breseq_" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json"
    output:
        dir = directory("breseq_" + sample_info.get_reference_prefix() + "/cov/{sample}"),
        #whole_genome_file = "evaluate/coverage_plots/{technology}/{sample}/all.pdf",
        done_file = "breseq_" + sample_info.get_reference_prefix() + "/cov/{sample}/coverage.done"
    log:
        "logs/coverage-plots-breseq_" + sample_info.get_reference_prefix() + "_{sample}.log"
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
