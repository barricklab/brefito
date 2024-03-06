try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "filter-nanopore-reads.smk"

rule all:
    input: [ "assemblies/{}.fasta".format(d) for d in sample_info.get_sample_list()]

rule assemble_with_flye:
    input:
        "nanopore-reads-filtered/{sample}.fastq"
    output:
        fasta = "assemblies/{sample}.fasta",
        gfa = "assemblies/{sample}.gfa",
        tmpdir = temp(directory("flye_assembly/{sample}"))
    log:
        "logs/flye_{sample}.log"
    conda:
        "../envs/flye.yml"
    threads: 8
    shell:
        """
        mkdir -p {output.tmpdir}
        flye --nano-raw {input} --threads {threads} --out-dir {output.tmpdir} > {log} 2>&1
        cp {output.tmpdir}/assembly.fasta {output.fasta}
        cp {output.tmpdir}/assembly_graph.gfa {output.gfa}
        """
