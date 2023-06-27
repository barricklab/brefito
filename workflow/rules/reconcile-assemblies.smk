import glob
import os.path


input_files=glob.glob("05_trycycler/*/cluster_*")
print(input_files)
OUTPUTS=[]
print("Running on the following datasets...")
for this_input_file in input_files:
    OUTPUTS.append(os.path.join(this_input_file, "2_all_seqs.fasta"))
print(OUTPUTS)

rule all:
    input: 
        OUTPUTS

rule trycycler_reconcile:
    input:
        reads = "01_trimmed_reads/{dataset}.long_reads.trimmed.fastq",
        cluster_dir = "05_trycycler/{dataset}/cluster_{cluster_id}"
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/reconcile.log"
    threads: 24
    shell:
        "trycycler reconcile --threads {threads} --reads {input.reads} --cluster_dir {input.cluster_dir} 2> {log}"

