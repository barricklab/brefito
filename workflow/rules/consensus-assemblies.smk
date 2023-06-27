import glob
import os.path


input_files=glob.glob("05_trycycler/*")
print(input_files)
OUTPUTS=[]
print("Running on the following datasets...")
for this_input_file in input_files:
    sample_name = os.path.basename(this_input_file) 
    OUTPUTS.append(os.path.join("06_trycycler_output", sample_name + ".fasta"))
print(OUTPUTS)

CLUSTER_IDS = ["001", "002", "003", "004"]

rule all:
    input: 
        OUTPUTS

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
        "06_trycycler_output/{dataset}.fasta"
    shell:
        "cat {input} > {output}"


rule trycycler_msa:
    input:
        all_seqs = "05_trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta",
        cluster_dir = "05_trycycler/{dataset}/cluster_{cluster_id}"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/msa.log"
    threads: 16
    shell:
        "trycycler msa --threads {threads} --cluster_dir {input.cluster_dir} 2> {log}"

rule trycycler_partition:
    input:
        reads = "01_trimmed_reads/{dataset}.long_reads.trimmed.fastq",
        all_seqs = "05_trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta",
        cluster_dir = "05_trycycler/{dataset}/cluster_{cluster_id}"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/partition.log"
    threads: 32
    shell:
        "trycycler partition --threads {threads} --reads {input.reads} --cluster_dir {input.cluster_dir} 2> {log}"

rule trycycler_consensus:
    input:
        all_seqs = "05_trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq",
        cluster_dir = "05_trycycler/{dataset}/cluster_{cluster_id}"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/7_final_consensus.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/consensus.log"
    threads: 32
    shell:
        "trycycler consensus --threads {threads} --cluster_dir {input.cluster_dir} 2> {log}"
