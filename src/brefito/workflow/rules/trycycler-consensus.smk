import glob
import os.path

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_trycycler_consensus:
    input:
        ["assemblies/" + s + ".fasta" for s in sample_info.get_sample_list() ]
    default_target: True

def find_all_contigs_for_sample(wildcards):
   from glob import glob
   path = "trycycler/{dataset}/cluster_{cluster_id}"
   files = glob(path.format(dataset=wildcards.dataset, cluster_id="*"))
   files = [s + "/7_final_consensus.fasta" for s in files]
   return files

rule cat_contigs:
    input:
        find_all_contigs_for_sample
    output:
        "assemblies/{dataset}.fasta"
    shell:
        """
        ARG={input}
        if [[ -z "$ARG" ]]; then
            echo "No inut files found for generating {output}"
            touch {output}
        else
            cat {input} > {output}
        fi
        
        """
rule trycycler_msa:
    input:
        "trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta",
    output:
        "trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "trycycler/{dataset}/cluster_{cluster_id}/msa.log"
    threads: 16
    shell:
        "trycycler msa --threads {threads} --cluster_dir trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

rule trycycler_partition:
    input:
        reads = "nanopore-reads-filtered/{dataset}.fastq",
        all_seqs = "trycycler/{dataset}/cluster_{cluster_id}/3_msa.fasta"
    output:
        "trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq"
    conda:
        "../envs/trycycler.yml"
    log: 
        "trycycler/{dataset}/cluster_{cluster_id}/partition.log"
    threads: 32
    shell:
        "trycycler partition --threads {threads} --reads {input.reads} --cluster_dir trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

rule trycycler_consensus:
    input:
        "trycycler/{dataset}/cluster_{cluster_id}/4_reads.fastq"
    output:
        "trycycler/{dataset}/cluster_{cluster_id}/7_final_consensus.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "trycycler/{dataset}/cluster_{cluster_id}/consensus.log"
    threads: 32
    shell:
        "trycycler consensus --threads {threads} --cluster_dir trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"
