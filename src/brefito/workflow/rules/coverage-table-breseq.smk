include: "predict-mutations-breseq.smk"

rule all_coverage_plots_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/tbl/" + s + ".done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to create  output files for each reference by default if no argument provided...
rule coverage_table_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        breseq_summary_json = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json"
    output:
        file = "breseq-" + sample_info.get_reference_prefix() + "/tbl/{sample}/",
        #whole_genome_file = "evaluate/coverage_plots/{technology}/{sample}/all.pdf",
        done_file = "breseq-" + sample_info.get_reference_prefix() + "/tbl/{sample}.done"
    log:
        "logs/coverage-table-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq bam2cov -b {input.bam} -f {input.fasta} -o {output.file} -t --resolution 0 -r REL606:1-4629812 -j {input.breseq_summary_json}  > {log} 2>&1
        touch {output.done_file}
        """
