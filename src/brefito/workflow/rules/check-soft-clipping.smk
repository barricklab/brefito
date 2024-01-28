include: "align-reads.smk"

BREFITO_PACKAGE_PATH = config["brefito_package_path"]

... ## Fix output file management
def all_soft_clipping_outputs():
    l = [] 
    for s in sample_info.get_sample_list():
        for r in sample_info.get_nanopore_read_base_list(s):
            l = l + ["aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/" + s + ".nanopore_reads." + r + ".summary.csv"]
        for r in sample_info.get_illumina_SE_read_base_list(s):
            l = l + ["aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/" + s + ".illumina_reads." + r + ".SE.summary.csv"]
        for r in sample_info.get_illumina_PE_read_base_list(s):
            l = l + ["aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/" + s + ".illumina_reads." + r + ".PE.summary.csv"]
    return(l)

rule all_evaluate_soft_clipping:
    input: all_soft_clipping_outputs()
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule tabulate_all_soft_clipping:
    input:
        bam = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam",
        fasta = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta",
    wildcard_constraints:
        technology = "nanopore_reads.*|illumina_reads.*.SE|illumina_reads.*.PE"
    output:
        csv = "aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/{sample}.{technology}.all.csv" 
    log:
        "logs/" + "check_soft_clipping_" + sample_info.get_reference_prefix() + "_{sample}.{technology}.log"
#    conda:
#        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq SOFT-CLIPPING -o {output.csv} -b {input.bam} -f {input.fasta} > {log} 2>&1
        """

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule summarize_soft_clipping:
    input:
        csv = "aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/{sample_technology}.all.csv" 
    output:
        plot = "aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/{sample_technology}.plot.pdf",
        summary = "aligned_reads_" + sample_info.get_reference_prefix() + "/soft_clipping/{sample_technology}.summary.csv" 
    log:
        "logs/" + "check_soft_clipping_summarize_" + sample_info.get_reference_prefix() + "_{sample_technology}.log"
    conda:
        "../envs/rstats.yml"
    threads: 1
    shell:
        """
        Rscript {BREFITO_PACKAGE_PATH}/summarize_soft_clipping.R -i {input.csv} -p {output.plot} -o {output.summary} > {log} 2>&1
        """

