include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

rule align_nanopore_reads:
    input:
        reads = "01_trimmed_nanopore_reads/{sample}.fastq",
        reference = "assemblies/{sample}.fasta"
    output:
        temp("intermediate_aligned_reads/{sample}/nanopore.sam")
    conda:
        "../envs/minimap2.yml"
    threads: 12
    shell:
        "minimap2 -t {threads} -x map-ont -a -o {output} {input.reference} {input.reads}"

rule align_illumina_reads:
    input:
        reads = expand("01_trimmed_illumina_reads/{{sample}}.R{read_num}.fastq.gz", read_num=["1", "2"]),
        reference = "assemblies/{sample}.fasta"
    output:
        temp("intermediate_aligned_reads/{sample}/illumina.sam")
    conda:
        "../envs/bowtie2.yml"
    threads: 12
    shell:
        """
        mkdir -p intermediate_aligned_reads/{wildcards.sample}
        bowtie2-build {input.reference} intermediate_aligned_reads/{wildcards.sample}/bowtie2
        bowtie2 --threads {threads} -x intermediate_aligned_reads/{wildcards.sample}/bowtie2 -1 {input.reads[0]} -2 {input.reads[1]} -S {output}
        rm intermediate_aligned_reads/{wildcards.sample}/bowtie2*
        """

rule samtools_view:
    input:
        "intermediate_aligned_reads/{sample}/{technology}.sam"
    output:
        temp("intermediate_aligned_reads/{sample}/{technology}.unsorted.bam")
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -bS {input} -o {output}
        """

rule samtools_sort:
    input:
        "intermediate_aligned_reads/{sample}/{technology}.unsorted.bam"
    output:
        "output_aligned_reads/{sample}/{technology}.bam"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools sort --threads {threads} {input} -o {output}
        """


rule samtools_bam_index:
    input:
        "output_aligned_reads/{sample}/{technology}.bam"
    output:
        "output_aligned_reads/{sample}/{technology}.bam.bai"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools index {input}
        """

rule samtools_fasta_index:
    input:
        "assemblies/{sample}.fasta"
    output:
        fasta = "output_aligned_reads/{sample}/reference.fasta",
        fai = "output_aligned_reads/{sample}/reference.fasta.fai"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        cp {input} {output.fasta}
        samtools faidx {output.fasta}
        """