# Does all steps through assembly and clustering

include: "load-data-csv.smk"

include: "filter-nanopore-reads.smk"
include: "subsample-illumina-reads.smk"


def get_short_reads(wildcards):
    sample = wildcards.sample
    illumina_read_list = []
    if (sample in illumina_read_lists.keys()):
        illumina_read_list = sorted(illumina_read_lists[sample])
        illumina_read_list = list(illumina_read_list[:2])
    return [ os.path.join("01_trimmed_illumina_reads", os.path.basename(s)) for s in illumina_read_list]

def get_long_reads(wildcards):
    sample = wildcards.sample
    nanopore_read_list = []
    if (sample in nanopore_read_lists.keys()):
        nanopore_read_list = nanopore_read_lists[sample]

    return [ os.path.join("01_trimmed_nanopore_reads", os.path.basename(s)) for s in nanopore_read_list]

def all_samples():
    return(['assemblies/' + item + ".fasta" for item in sample_list])

rule all:
    input: all_samples()

rule assemble_with_unicycler:
    input:
        short_reads = get_short_reads,
        long_reads = get_long_reads
    output:
        fasta = "assemblies/{sample}.fasta",
        gfa = "assemblies/{sample}.gfa",
        directory = temp(directory("03_unicycler_assembly_temp/{sample}"))
    log:
        "logs/unicycler_{sample}.log"
    conda:
        "../envs/unicycler.yml"
    threads: 4
    shell:
        """
        if [ -n "{input.long_reads}" ]; then
          unicycler --threads {threads} -1 {input.short_reads[0]} -2 {input.short_reads[1]} -l {input.long_reads} -o {output.directory} > {log} 2>&1
        else
          unicycler --threads {threads} -1 {input.short_reads[0]} -2 {input.short_reads[1]} -o {output.directory} > {log} 2>&1
        fi
        cp {output.directory}/assembly.fasta {output.fasta}
        cp {output.directory}/assembly.gfa {output.gfa}
        """