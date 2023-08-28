import glob
import os.path

def find_available_files(wildcards):
   from glob import glob
   path = "05_trycycler/{dataset}/cluster_{cluster_id}"
   files = glob(path.format(dataset=wildcards.dataset, cluster_id="*"))
   files = [s + "/7_final_consensus.fasta" for s in files]
   return files

rule cat_contigs:
    input:
        find_available_files
    output:
        "assemblies/{dataset}.fasta"
    shell:
        "cat {input} > {output}"

rule trycycler_msa:
    input:
        "05_trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta",
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/msa.log"
    threads: 16
    shell:
        "trycycler msa --threads {threads} --cluster_dir 05_trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

rule trycycler_partition:
    input:
        reads = "02_filtered_nanopore_reads/{dataset}.fastq",
        all_seqs = "05_trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/partition.log"
    threads: 32
    shell:
        "trycycler partition --threads {threads} --reads {input.reads} --cluster_dir 05_trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

rule trycycler_consensus:
    input:
        "05_trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/7_final_consensus.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/consensus.log"
    threads: 32
    shell:
        "trycycler consensus --threads {threads} --cluster_dir 05_trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"
