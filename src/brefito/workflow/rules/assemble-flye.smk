include: "filter-nanopore-reads.smk"

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


rule assemble_with_flye:
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_01.fastq"
    output:
        fasta = "assemblies/{dataset}.fasta",
        gfa = "assemblies/{dataset}.gfa",
    log:
        "logs/flye_{dataset}.log"
    params:
        tmpdir = "04_flye_assembly/{dataset}"
    conda:
        "../envs/flye.yml"
    threads: 16
    shell:
        """
        mkdir -p {params.tmpdir}
        flye --nano-raw {input} --threads {threads} --out-dir {params.tmpdir} 2> {log}
        cp {params.tmpdir}/assembly.fasta {output.fasta}
        cp {params.tmpdir}/assembly_graph.gfa {output.gfa}
        """
