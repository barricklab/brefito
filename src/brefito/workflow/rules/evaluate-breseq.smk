include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

BRESEQ_OPTIONS = ""
if "breseq_options" in config.keys():
    BRESEQ_OPTIONS = config["breseq_options"]

rule evaluate_breseq_with_nanopore_reads:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = "01_trimmed_nanopore_reads/{sample}.fastq.gz"
    output:
        intermediate_path = directory("06_breseq/{sample}_nanopore_reads"),
        output_path = directory("evaluate/breseq_output/{sample}_nanopore_reads"),
        data_path = directory("evaluate/breseq_data/{sample}_nanopore_reads")
    log: 
        "logs/evaluate_breseq_nanopore_reads_{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 8
    shell:
        """
        breseq --nanopore {BRESEQ_OPTIONS} -j {threads} -r {input.reference} -o {output.intermediate_path} {input.reads} > {log} 2>&1
        mv {output.intermediate_path}/output {output.output_path}
        mv {output.intermediate_path}/data {output.data_path}
       """    

rule evaluate_breseq_with_illumina_reads:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = expand("01_trimmed_illumina_reads/{{sample}}.R{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        intermediate_path = directory("06_breseq/{sample}_illumina_reads"),
        output_path = directory("evaluate/breseq_output/{sample}_illumina_reads"),
        data_path = directory("evaluate/breseq_data/{sample}_nanopore_reads")
    log: 
        "logs/evaluate_breseq_illumina_reads_{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 8
    shell:
        """
        breseq -j {threads} {BRESEQ_OPTIONS} -r {input.reference} -o {output.intermediate_path}  {input.reads} > {log} 2>&1
        mv {output.intermediate_path}/output {output.output_path}
        mv {output.intermediate_path}/data {output.data_path}
        """    