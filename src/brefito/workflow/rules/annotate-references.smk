rule annotate_with_prokka:
    input:
        "references/{sample}.fasta"
    output:
        directory = directory("01_prokka/{sample}")
        annotated_reference = "01_prokka/{sample}/reference.gbk"
    conda:
        "../envs/prokka.yml"
    threads: 12
    shell:
        "prokka --cpus {threads} --outdir {output.directory} {input} --count 2 --out_dir 03_subsampled_nanopore_reads/{wildcards.sample}"


rule annotate_with_isescan:
    input:
        "references/{sample}.fasta"
    output:
        fasta = "assemblies/{dataset}.fasta",
        gfa = "assemblies/{dataset}.gfa",
    log:
        "logs/flye_{dataset}.log"
    params:
        tmpdir = "04_flye_assembly/{dataset}"
    conda:
        "../envs/isescan.yml"
    threads: 16
    shell:
        """
        mkdir -p {params.tmpdir}
        flye --nano-raw {input} --threads {threads} --out-dir {params.tmpdir} 2> {log}
        cp {params.tmpdir}/assembly.fasta {output.fasta}
        cp {params.tmpdir}/assembly_graph.gfa {output.gfa}
        """

rule combine_annotation_with_breseq:
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
        "../envs/breseq.yml"
    threads: 16
    shell:
        """
        breseq CONVERT-REFEENCE
        """
