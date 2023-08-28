BRESEQ_OPTIONS = ""
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

rule trycycler_reconcile:
    input:
        reads = "02_filtered_nanopore_reads/{dataset}.fastq",
    output:
        "05_trycycler/{dataset}/cluster_{cluster_id}/2_all_seqs.fasta"
    conda:
        "../envs/trycycler.yml"
    log: 
        "05_trycycler/{dataset}/cluster_{cluster_id}/reconcile.log"
    threads: 12
    shell:
        "trycycler reconcile {TRYCYCLER_OPTIONS} --threads {threads} --reads {input.reads} --cluster_dir 05_trycycler/{wildcards.dataset}/cluster_{wildcards.cluster_id} 2> {log}"

