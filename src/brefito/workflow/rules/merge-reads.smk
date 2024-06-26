try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "download-data.smk"

rule merge_all_reads:
    input:
        nanopore_reads = ["merged-reads/{}.nanopore.fastq.gz".format(s) for s in sample_info.get_samples_with_nanopore_reads()],
        illumina_SE_reads = ["merged-reads/{}.illumina.SE.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_SE_reads()],
        illumina_R1_reads = ["merged-reads/{}.illumina.R1.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_PE_reads()],
        illumina_R2_reads = ["merged-reads/{}.illumina.R2.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_PE_reads()]
    default_target: True

rule merge_nanopore_reads:
    input:
        lambda wildcards: [ "nanopore-reads/" + d for d in sample_info.get_nanopore_read_list(wildcards.sample)]
    output:
        "merged-reads/{sample}.nanopore.fastq.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule merge_illumina_reads:
    input:
        lambda wildcards: [ "illumina-reads/" + d for d in sample_info.get_file_list(wildcards.sample, "illumina-" + wildcards.illumina_type)]
    output:
        "merged-reads/{sample}.illumina.{illumina_type}.fastq.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """
