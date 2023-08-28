include: "trim-nanopore-reads.smk"

rule polish_with_medaka:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = "02_filtered_nanopore_reads/{sample}.fastq"
    output:
        path = temp(directory("06_medaka_polish/{sample}")),
        fasta = "assemblies/{sample}.fasta.polished"
    log: 
        "logs/polish_medaka_{sample}.log"
    conda:
        "../envs/medaka.yml"
    threads: 8
    shell:
      """
      # Copy fasta to output directory in order to avoid preserving an index file
      mkdir -p {output.path}
      cp {input.reference} {output.path}
      medaka_consensus -i {input.reads} -d {output.path}/{wildcards.sample}.fasta -o {output.path} -m r941_min_sup_g507 -t {threads} > {log} 2>&1
      mv {output.path}/consensus.fasta {output.fasta}
      """    