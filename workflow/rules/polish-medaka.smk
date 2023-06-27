import glob
DATASETS=glob.glob("input/*.fasta")
DATASETS = [s.replace('input/', 'output/') for s in DATASETS]
print(DATASETS)

rule all:
  input:
        DATASETS

rule polish_with_medaka:
    input:
        reference = "input/{sample}.fasta",
        reads = "01_trimmed_reads/{sample}.long_reads.trimmed.fastq"
    output:
        path = directory("06_medaka_polish/{sample}"),
        fasta = "output/{sample}.fasta"
    log: 
        "06_medaka_polish/{sample}.log"
    conda:
        "../envs/medaka.yml"
    threads: 8
    shell:
      """
      medaka_consensus -i {input.reads} -d {input.reference} -o {output.path} -m r941_min_sup_g507 -t {threads} 2> {log}
      mv {output.path}/consensus.fasta {output.fasta}
      """    