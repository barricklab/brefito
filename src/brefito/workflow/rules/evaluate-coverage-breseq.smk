try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "predict-mutations-breseq.smk"


rule all_evaluate_coverage:
    input:
        ["breseq-references/coverage_plots/" + s + "/coverage.done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule evaluate_coverage:
    input:
        bam = "breseq-references/data/{sample}/data/reference.bam",
        bai = "breseq-references/data/{sample}/data/reference.bam.bai",
        fasta = "breseq-references/data/{sample}/data/reference.fasta",
        fai = "breseq-references/data/{sample}/data/reference.fasta.fai",
        breseq_summary_json = "breseq-references/data/{sample}/data/summary.json"
    output:
#        dir = directory("breseq-references/coverage_plots/{sample}"),
        done_file = "breseq-references/coverage_plots/{sample}/coverage.done" 
    log:
        "logs/evaluate_coverage_{sample}.log"
#    Official breseq doesn't support plotting entire sequence (yet)
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        mkdir -p breseq-references/coverage_plots/{wildcards.sample}
        # tiles across sequences
        breseq bam2cov --format PNG -o {output.dir} -b {input.bam} -f {input.fasta} -a --tile-size 10000 --tile-overlap 2000 -p 1000 -j {input.breseq_summary_json} > {log} 2>&1
        # entire sequences
        breseq bam2cov --format PNG -o {output.dir}/{wildcards.sample} -b {input.bam} -f {input.fasta} -a -p 10000 -j {input.breseq_summary_json} > {log} 2>&1
        touch {output.done_file}
        """
        #breseq bam2cov -o {output.whole_genome_file} -b {input.bam} -f {input.fasta} -a -j {input.breseq_summary_json} -p 10000 > {log} 2>&1
