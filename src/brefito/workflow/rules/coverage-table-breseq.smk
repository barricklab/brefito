include: "predict-mutations-breseq.smk"

rule all_coverage_plots_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/coverage_table/" + s + "/coverage_table.done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to create  output files for each reference by default if no argument provided...
rule coverage_table_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
    output:
        output_path = directory("breseq-" + sample_info.get_reference_prefix() + "/coverage_table/{sample}"),
        done_file = "breseq-" + sample_info.get_reference_prefix() + "/coverage_table/{sample}/coverage_table.done"
    log:
        "logs/coverage-table-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq bam2cov -b {input.bam} -f {input.fasta} -o {output.output_path} -t --resolution 0 > {log} 2>&1
        touch {output.done_file}
        """
