try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"
include: "merge-references.smk"

def get_all_output_files_list():
    all_output_files_list = []

    for s in sample_info.get_sample_list():

        all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.fasta")
        all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.fasta.fai")
        all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.gff3")
        
        for n in sample_info.get_nanopore_read_base_list(s):
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads.{}.bam".format(n))
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads.{}.bam.bai".format(n))

        for i_se in sample_info.get_illumina_SE_read_base_list(s):
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.SE.bam".format(i_se))
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.SE.bam.bai".format(i_se))

        for i_pe in sample_info.get_illumina_PE_read_base_list(s):
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.PE.bam".format(i_pe))
            all_output_files_list.append("aligned_reads_" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.PE.bam.bai".format(i_pe))

    print(all_output_files_list)

    return all_output_files_list

rule all_align_reads:
    input:
        get_all_output_files_list()
    default_target: True

rule align_nanopore_reads:
    input:
        reads = "nanopore_reads_trimmed/{reads}.fastq.gz",
        reference = sample_info.get_merged_reference_prefix() + "/{sample}.fasta"
    output:
        temp("aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/nanopore_reads.{reads}.sam")
    conda:
        "../envs/minimap2.yml"
    threads: 12
    shell:
        "minimap2 --secondary=no -t {threads} -x map-ont -a -o {output} {input.reference} {input.reads}"

rule align_PE_illumina_reads:
    input:
        reads = expand("illumina_reads_trimmed/{{reads}}.{read_num}.fastq.gz", read_num=["R1", "R2"]),
        reference = sample_info.get_merged_reference_prefix() + "/{sample}.fasta"
    output:
        temp("aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.PE.sam")
    conda:
        "../envs/bowtie2.yml"
    params:
        bowtie2_index = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.PE.bowtie2_index"
    threads: 12
    shell:
        """
        mkdir -p aligned_reads/data/{wildcards.sample}
        bowtie2-build {input.reference} {params.bowtie2_index}
        bowtie2 --threads {threads} -x {params.bowtie2_index} -1 {input.reads[0]} -2 {input.reads[1]} -S {output}
        rm {params.bowtie2_index}*
        """

rule align_SE_illumina_reads:
    input:
        reads = expand("illumina_reads_trimmed/{{reads}}.{read_num}.fastq.gz", read_num=["SE"]),
        reference = sample_info.get_merged_reference_prefix() + "/{sample}.fasta"
    output:
        temp("aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.SE.sam")
    conda:
        "../envs/bowtie2.yml"
    params:
        bowtie2_index = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.SE.bowtie2_index"
    threads: 12
    shell:
        """
        mkdir -p aligned_reads/data/{wildcards.sample}
        bowtie2-build {input.reference} {params.bowtie2_index}
        bowtie2 --threads {threads} -x {params.bowtie2_index} -U {input.reads} -S {output}
        rm {params.bowtie2_index}*
        """

rule samtools_view:
    input:
        "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.sam"
    output:
        temp("aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.unsorted.bam")
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -bS {input} -o {output}
        """

rule samtools_sort:
    input:
        "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.unsorted.bam"
    output:
        "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools sort --threads {threads} {input} -o {output}
        """


rule samtools_bam_index:
    input:
        "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam"
    output:
        "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam.bai"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools index {input}
        """

rule copy_references:
    input:
        fasta = "merged_references/{sample}.fasta",
        gff3 = "merged_references/{sample}.gff3"
    output:
        fasta = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta",
        gff3 = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/reference.gff3"
    shell:
        """
        cp {input.fasta} {output.fasta}
        cp {input.gff3} {output.gff3}
        """

rule samtools_fasta_index:
    input:
        fasta = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta"
    output:
        fai = "aligned_reads_" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta.fai"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools faidx {input.fasta}
        """