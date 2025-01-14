#!/usr/bin/env bash
# Author: Ira Zibbu
# A bash script that accepts the name of a sample, and returns the consensus assembly
# This script performs clustering/ trimming/ dot plots/ resolving/ combining
# Usage: autocycler_cluster_trim_resolve_combine.sh <name of sample>

sample=$1
input_directory="autocycler/$sample"

# run cluster
autocycler cluster -a ${input_directory}

for c in "${input_directory}/clustering/qc_pass/cluster_"*; do
    autocycler trim -c "$c"
    autocycler dotplot -i "$c"/1_untrimmed.gfa -o "$c"/1_untrimmed.png
    autocycler dotplot -i "$c"/2_trimmed.gfa -o "$c"/2_trimmed.png
    autocycler resolve -c "$c"
done

output_directory="autocycler/$sample/output"

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi

autocycler combine -a "${output_directory}" -i ${input_directory}/clustering/qc_pass/cluster_*/5_final.gfa
