# Commands for testing...

# Clean up prior results
rm -rf breseq-* 
rm -rf illumina-* 
rm -rf nanopore-* g
rm -rf genome_diffs 
rm -rf references 
rm -rf mutants
rm -rf aligned-reads-* 
rm -rf logs

# Run various commands on references specified in data.csv
brefito predict-mutations-breseq --config breseq_options="-l 50" --resources "connections=4"
brefito coverage-plots-breseq
brefito align-reads
brefito check-soft-clipping

# Run breseq on the genome with the predicted mutations applied
#Copies over the generated genome_diff so we can apply it
mkdir -p genome_diffs
cp breseq_references/gd/Ara-1_50000gen_11331.gd genome_diffs
brefito mutate-genomes-gdtools

# Run commands on new reference generated in "mutants"
brefito predict-mutations-breseq-mutants --config breseq_options="-l 50"
brefito coverage-plots-breseq-mutants
brefito align-reads-mutants
brefito check-soft-clipping-mutants

# Run commands on reference in "assemblies"
brefito predict-mutations-breseq-assemblies --config breseq_options="-l 50"
brefito coverage-plots-breseq-assemblies
brefito align-reads-assemblies
brefito check-soft-clipping-assemblies

