try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-illumina-reads.smk"


rule all_polish_polypolish:
    input:
        ["assemblies/" + s + ".fasta.polished" for s in sample_info.get_sample_list() ]
    default_target: True

# Merge the trimmed PE reads so that we use all that were provided
rule merge_PE_illumina_reads:
    input:
        lambda wildcards: [ "illumina-reads-trimmed/" + d for d in sample_info.get_file_list(wildcards.sample, "illumina-R" + wildcards.read_num)]
    output:
        "reads-trimmed-merged/{sample}.illumina.R{read_num}.fastq.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """


rule polish_with_polypolish:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = expand("reads-trimmed-merged/{{sample}}.illumina.R{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        path = temp(directory("polypolish/{sample}")),
        polished_fasta = "assemblies/{sample}.fasta.polished"
    conda:
        "../envs/polypolish.yml"
    log:
        "logs/polish_polypolish_{sample}.log"
    threads: 16 
    shell:
      """
      mkdir -p {output.path}
      cp {input.reference} {output.path}/reference.fasta 
      bwa index {output.path}/reference.fasta 2> {log}
      bwa mem -t 16 -a {output.path}/reference.fasta {input.reads[0]} > {output.path}/alignments_1.sam 2>> {log}
      bwa mem -t 16 -a {output.path}/reference.fasta {input.reads[1]} > {output.path}/alignments_2.sam 2>> {log}
      polypolish polish {output.path}/reference.fasta {output.path}/alignments_1.sam {output.path}/alignments_2.sam > {output.polished_fasta} 2>> {log}
      """    
