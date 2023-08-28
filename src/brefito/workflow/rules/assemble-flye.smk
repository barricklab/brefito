include: "filter-nanopore-reads.smk"

def get_fastq_reads(wildcards):
    #print(wildcards)
    if ("subsample" in config.keys()) and (config[subsample] != False):
        return ["03_subsampled_nanopore_reads/{}/sample_01.fastq".format(wildcards.dataset)]
    else:
        return ["02_filtered_nanopore_reads/{}.fastq".format(wildcards.dataset)]
    return []

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
        get_fastq_reads
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
