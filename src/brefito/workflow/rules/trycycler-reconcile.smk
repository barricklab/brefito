import glob
import os.path

if "trycycler_options" in config.keys():
    TRYCYCLER_OPTIONS = config["trycycler_options"]
else:
    TRYCYCLER_OPTIONS = "--min_1kbp_identity 15 --max_add_seq 50000 --max_add_seq_percent 5"

print ()
print ("Using these command-line options for trycycler reconcile:")
print ("   " + TRYCYCLER_OPTIONS)
print ("Add something like this to your brefito command line to override these options:")
print ("   --config \"trycycler_options= --min_1kbp_identity 15 --max_add_seq 50000 --max_add_seq_percent 5\"")
print()

def get_all_trycycler_reconcile_outputs():
    smk_targets = []
    input_files=glob.glob("trycycler/*/cluster_*")
    for this_input_file in input_files:
        smk_targets.append(os.path.join(this_input_file, "2_all_seqs.fasta"))
    return(smk_targets)

rule all_trycycler_reconcile:
    input:
        get_all_trycycler_reconcile_outputs()
    default_target: True

rule trycycler_reconcile:
    input:
        reads = "nanopore-reads-filtered/{dataset}.fastq",
    output:
        "trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "trycycler/{dataset}/cluster_{cluster_id}/reconcile.log"
    threads: 12
    shell:
        "trycycler reconcile {TRYCYCLER_OPTIONS} --threads {threads} --reads {input.reads} --cluster_dir trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

