try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

rule merge_all_trimmed_reads:
    input:
        # long reads of every type: <sample>.<type>, e.g. <sample>.nanopore or <sample>.pacbio-hifi
        long_reads = ["merged-reads-trimmed/{}.{}.fastq.gz".format(s, sample_info.get_long_read_type(s)) for s in sample_info.get_samples_with_long_reads()],
        illumina_SE_reads = ["merged_reads_trimmed/{}.illumina.SE.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_SE_reads()],
        illumina_R1_reads = ["merged_reads_trimmed/{}.illumina.R1.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_PE_reads()],
        illumina_R2_reads = ["merged_reads_trimmed/{}.illumina.R2.fastq.gz".format(s) for s in sample_info.get_samples_with_illumina_PE_reads()]
    default_target: True

rule merge_nanopore_reads:
    input:
        lambda wildcards: [ wildcards.lrt + "-reads-trimmed/" + d for d in sample_info.get_long_read_list(wildcards.sample)]
    output:
        "merged-reads-trimmed/{sample}.{lrt}.fastq.gz"
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule merge_illumina_reads:
    input:
        lambda wildcards: [ "illumina_reads_trimmed/" + d for d in sample_info.get_file_list(wildcards.sample, "illumina-" + wildcards.illumina_type)]
    output:
        "merged_reads_trimmed/{sample}.illumina.{illumina_type}.fastq.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """
