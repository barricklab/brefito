include: "align-reads.smk"

BREFITO_PACKAGE_PATH = brefito_config["BREFITO_PACKAGE_PATH"]

rule all_evaluate_soft_clipping:
    input: 
        ["aligned-reads-" + sample_info.get_reference_prefix() + "/soft-clipping/" + s + ".plot.pdf" for s in sample_info.get_sample_list()]
    default_target: True


def get_bam_for_sample(sample):
    bam_list = []
    for r in sample_info.get_nanopore_read_base_list(sample):
        bam_list += ["aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + sample + "/nanopore_reads." + r + ".bam"]
    for r in sample_info.get_illumina_SE_read_base_list(sample):
        bam_list += ["aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + sample + "/illumina_reads." + r + ".SE.bam"]
    for r in sample_info.get_illumina_PE_read_base_list(sample):
        bam_list += ["aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + sample + "/illumina_reads." + r + "PE.bam"]
    return (bam_list)

# We use the breseq files to get the average coverage
# Need to update bam2cov to by default output files for each reference...
rule tabulate_all_soft_clipping:
    input:
        bam = lambda wildcards: get_bam_for_sample(wildcards.sample),
        fasta = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta",
        fai = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta.fai",
    output:
        csv = "aligned-reads-" + sample_info.get_reference_prefix() + "/soft-clipping/{sample}.all.csv" 
    log:
        "logs/" + "check-soft-clipping-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq SOFT-CLIPPING -o {output.csv} -f {input.fasta} {input.bam} > {log} 2>&1
        """

rule summarize_soft_clipping:
    input:
        csv = "aligned-reads-" + sample_info.get_reference_prefix() + "/soft-clipping/{sample}.all.csv" 
    output:
        plot = "aligned-reads-" + sample_info.get_reference_prefix() + "/soft-clipping/{sample}.plot.pdf",
        summary = "aligned-reads-" + sample_info.get_reference_prefix() + "/soft-clipping/{sample}.summary.csv" 
    log:
        "logs/" + "check-soft-clipping-summarize-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/rstats.yml"
    threads: 1
    shell:
        """
        Rscript {BREFITO_PACKAGE_PATH}/summarize_soft_clipping.R -i {input.csv} -p {output.plot} -o {output.summary} > {log} 2>&1
        """

