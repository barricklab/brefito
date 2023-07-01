
rule trycycler_reconcile:
    input:
        reads = "02_filtered_nanopore_reads/{dataset}.fastq",
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

