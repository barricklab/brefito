try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "filter-nanopore-reads.smk"


if 'FLYE_OPTIONS' in brefito_config.keys():
    FLYE_OPTIONS = brefito_config['FLYE_OPTIONS']
    print("flye assembly with user-specified command-line options:")
    print("  " + FLYE_OPTIONS)
else:
    FLYE_OPTIONS="--nano-raw"
    print("flye assembly with default command-line options:")
    print("  " + FLYE_OPTIONS)
    print("To change these options add something like this to your brefito command")
    print("  --config FLYE_OPTIONS=\"--meta --nano-hq")

rule all:
    input: [ "assemblies/{}.fasta".format(d) for d in sample_info.get_samples_with_nanopore_reads()]

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
        flye {FLYE_OPTIONS} {input} --threads {threads} --out-dir {output.tmpdir} > {log} 2>&1
        cp {output.tmpdir}/assembly.fasta {output.fasta}
        cp {output.tmpdir}/assembly_graph.gfa {output.gfa}
        """
