# Does all steps through assembly and clustering

include: "filter-nanopore-reads.smk"
include: "subsample-illumina-reads.smk"

READ_NUMS = ["1", "2"]

rule subsample_nanopore_reads:
    input:
        "02_filtered_nanopore_reads/{sample}.fastq"
    output:
        "03_subsampled_nanopore_reads/{sample}/sample_01.fastq",
        "03_subsampled_nanopore_reads/{sample}/sample_02.fastq"
    conda:
        "../envs/trycycler.yml"
    threads: 1
    shell:
        "trycycler subsample --genome_size {config[genome_size]} --reads {input} --count 2 --out_dir 03_subsampled_nanopore_reads/{wildcards.sample}"

rule assemble_with_unicycler:
    priority: 1
    input:
        long_reads = "03_subsampled_nanopore_reads/{dataset}/sample_01.fastq",
        short_reads = expand("03_subsampled_illumina_reads/{{dataset}}/sample_01_R{read_num}.fastq.gz", read_num=READ_NUMS)
    wildcard_constraints:
        assembly_id="01",
    output:
        fasta = "assemblies_unicycler/{dataset}.fasta",
        gfa = "assemblies_unicycler/{dataset}.gfa",
        directory = directory("03_unicycler_assembly_temp/{dataset}")
    conda:
        "../envs/unicycler.yml"
    threads: 24
    shell:
        """
        unicycler --threads {threads} -1 {input.short_reads[0]} -2 {input.short_reads[1]} -l {input.long_reads} -o {output.directory}
        cp {output.directory}/assembly.fasta {output.fasta}
        cp {output.directory}/assembly.gfa {output.gfa}
        """