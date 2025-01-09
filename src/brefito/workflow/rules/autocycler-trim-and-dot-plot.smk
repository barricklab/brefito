try: sample_info
except NameError:
    include: "load-sample-info.smk"

include: "autocycler-assemble.smk"

rule all_targets_autocycler_cluster_and_dotplot:
    input:
        ["autocycler/" + s + "/clustering/qc_pass/cluster_001/2_trimmed.gfa" for s in sample_info.get_sample_list()]
        ["autocycler/" + s + "/clustering/qc_pass/cluster_001/1_untrimmed.png" for s in sample_info.get_sample_list()]
        ["autocycler/" + s + "/clustering/qc_pass/cluster_001/2_trimmed.png" for s in sample_info.get_sample_list()]

    default_target: True

# For each sample, autocycler cluster creater >1 cluster, but the number of clusters will vary per sample
# Following target rules have been written such that snakemake only looks for the input-output dependencies for just one cluster
# But uses a bash loop to execute the trimming steps for all of the cluster for a sample

rule autocycler_trim_clusters:
    input:
        "autocycler/{dataset}/clustering/qc_pass/cluster_001/1_untrimmed.gfa"
    output:
        "autocycler/{dataset}/clustering/qc_pass/cluster_001/2_trimmed.gfa"
    threads: 16
    log:
        "logs/{dataset}/autocycler_trim_clusters.log"
    shell:
        """
        for c in autocycler/{wildcards.dataset}/clustering/qc_pass/cluster_*; do
        autocycler trim -c $c
        done
        """

# autocycler dotplot will make accept a gfa file for each cluster and make a dotplot for it
# for the same reasons as above, snakemake is only going to look for the input-output dependencies for the first cluster
# but use a bash shell for loop to run the commands for all of the clusters for a given sample

rule autocycler_make_dot_plots:
    input:
        trimmed_gfa = "autocycler/{dataset}/clustering/qc_pass/cluster_001/2_trimmed.gfa"
    output:
        untrimmed_dotplots = "autocycler/{dataset}/clustering/qc_pass/cluster_001/1_untrimmed.png"
        trimmed_dotplots = "autocycler/{dataset}/clustering/qc_pass/cluster_001/2_trimmed.png"
    params:
        untrimmed_gfa = "autocycler/{dataset}/clustering/qc_pass/cluster_001/1_untrimmed.gfa"
    log:
        "logs/{dataset}/autocycler_trim_clusters.log"
    shell:
        """
        for c in autocycler/{wildcards.dataset}/clustering/qc_pass/cluster_*; do
        cd $c
        autocycler dotplot -i 1_untrimmed.gfa -o 1_untrimmed.png
        autocycler gotplot -i 2_trimmed.gfa -o 2_trimmed.png
        cd ../
        done
        """
