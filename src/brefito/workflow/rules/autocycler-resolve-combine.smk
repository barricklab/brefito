try: sample_info
except NameError:
    include: "load-sample-info.smk"

include: "autocycler-trim-and-dot-plot.smk"

rule all_target_combine_assembly:
    input:
        ["autocycler/" + s + "/output/consensus_assembly.fasta" for s in sample_info.get_sample_list()]

# For each sample, autocycler cluster creater >=1 cluster, but the number of clusters will vary per sample
# Following target rules have been written such that snakemake only looks for the input-output dependencies for just one cluster
# But uses a bash loop to execute the trimming steps for all of the cluster for a sample

rule autocycler_resolve:
    input:
        "autocycler/{dataset}/clustering/qc_pass/cluster_001/2_trimmed.gfa"
    output:
        "autocycler/{dataset}/clustering/qc_pass/cluster_001/5_final.gfa"
    threads: 16
    log:
        "logs/{dataset}/autocycler_resolve.log"
    shell:
        """
        for c in autocycler/{wildcards.dataset}/clustering/qc_pass/cluster_*; do
        autocycler resolve -c $c
        done
        """

rule autocycler_combine:
    input:
        "autocycler/{dataset}/clustering/qc_pass/cluster_001/5_final.gfa"
    output:
        "autocycler/{dataset}/output/consensus_assembly.fasta"
    threads: 16
    log:
        "logs/{dataset}/autocycler_combine.log"
    shell:
        """
        autocycler combine -a autocycler/{wildcards.dataset}/output -i autocycler/{dataset}/clustering/qc_pass/cluster_*/5_final.gfa"
        """
