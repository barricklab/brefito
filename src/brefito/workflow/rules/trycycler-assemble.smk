# Does all steps through assembly and clustering

include: "subsample-nanopore-reads.smk"

num_resampled_read_sets=24
FLYE_RESAMPLED_IDS = ["01", "07", "13", "19", "25"]
MINIASM_MINIPOLISH_RESAMPLED_IDS = ["02", "08", "14", "20", "26"]
RAVEN_RESAMPLED_IDS = ["03", "09", "15", "21", "27"]
CANU_RESAMPLED_IDS = ["04", "10", "16", "22", "28"]
NECAT_RESAMPLED_IDS = ["05", "11", "17", "23", "29"]
UNICYCLER_RESAMPLED_IDS = []

FULL_FASTQ_IDS = ["25", "26", "27", "28", "29"]

# Uncomment if using unicycler
#include: "subsample-illumina-reads.smk"
#UNICYCLER_RESAMPLED_IDS = ["06", "12", "18", "24"]

ALL_RESAMPLED_IDS = FLYE_RESAMPLED_IDS + RAVEN_RESAMPLED_IDS + MINIASM_MINIPOLISH_RESAMPLED_IDS + CANU_RESAMPLED_IDS + NECAT_RESAMPLED_IDS + UNICYCLER_RESAMPLED_IDS

# Could copy but trying ln to save space
rule copy_full_fastq:
    input:
        "02_filtered_nanopore_reads/{dataset}.fastq"
    wildcard_constraints:
        assembly_id="|".join(FULL_FASTQ_IDS),
    output:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    threads: 1
    shell:
        "cp {input} {output}"


rule assemble_with_flye:
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(FLYE_RESAMPLED_IDS),
    output:
        "04_assemblies/{dataset}/assembly_{assembly_id}.fasta"
    log:
        "logs/{dataset}/flye_{assembly_id}.log"
    conda:
        "../envs/flye.yml"
    threads: 16
    shell:
        "mkdir -p 03_assembly_intermediates/{wildcards.dataset}/assembly_{wildcards.assembly_id} && flye --nano-raw {input} --threads {threads} --out-dir 03_assembly_intermediates/{wildcards.dataset}/assembly_{wildcards.assembly_id} > {log} 2>&1 && cp 03_assembly_intermediates/{wildcards.dataset}/assembly_{wildcards.assembly_id}/assembly.fasta {output} && cp 03_assembly_intermediates/{wildcards.dataset}/assembly_{wildcards.assembly_id}/assembly_graph.gfa 03_assembly_intermediates/{wildcards.dataset}/assembly_{wildcards.assembly_id}_01.gfa"


rule assemble_with_miniasm_minipolish:
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(MINIASM_MINIPOLISH_RESAMPLED_IDS),
    output:
        assembly = "04_assemblies/{dataset}/assembly_{assembly_id}.fasta",
        graph = "04_assemblies/{dataset}/assembly_{assembly_id}.gfa"
    log:
        "logs/{dataset}/miniasm_minipolish_{assembly_id}.log"
    conda:
        "../envs/miniasm_minipolish.yml"
    threads: 16
    shell:
        """
        # Create temporary intermediate files.
        overlaps=$(mktemp)".paf"
        unpolished_assembly=$(mktemp)".gfa"

        # Find read overlaps with minimap2.
        minimap2 -x ava-ont -t {threads} {input} {input} > "$overlaps" 2>> {log}

        # Run miniasm to make an unpolished assembly.
        miniasm -f {input} "$overlaps" > "$unpolished_assembly" 2>> {log}

        # Polish the assembly with minipolish, outputting the result to stdout.
        minipolish --threads {threads} {input} "$unpolished_assembly" > {output.graph} 2>> {log}

        # Convert to fasta
        any2fasta {output.graph} > {output.assembly} 2>> {log}
        
        # Clean up.
        rm "$overlaps" "$unpolished_assembly"
        """

rule assemble_with_raven:
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(RAVEN_RESAMPLED_IDS),
    output:
        "04_assemblies/{dataset}/assembly_{assembly_id}.fasta"
    log:
        "logs/{dataset}/raven_{assembly_id}.log"
    conda:
        "../envs/raven.yml"
    threads: 16
    shell:
        "raven --threads {threads} --disable-checkpoints --graphical-fragment-assembly 04_assemblies/{wildcards.dataset}/assembly_{wildcards.assembly_id}.gfa {input} > {output} 2> {log}"


rule assemble_with_canu:
    priority: 2
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(CANU_RESAMPLED_IDS),
    output:
        assembly_fasta = "03_canu_assembly_temp/{dataset}/{assembly_id}/canu.contigs.fasta",
        assembly_directory = directory("03_canu_assembly_temp/{dataset}/{assembly_id}")
    conda:
        "../envs/canu.yml"
    threads: 16
    shell:
        """
        canu -p canu -d {output.assembly_directory} -fast genomeSize="{config[genome_size]}" useGrid=false minThreads="{threads}" maxThreads="{threads}" -nanopore {input}
        rm -rf canu_temp
        """

rule trim_canu_assembly:
    priority: 3
    input:
        "03_canu_assembly_temp/{dataset}/{assembly_id}/canu.contigs.fasta"
    wildcard_constraints:
        assembly_id="|".join(CANU_RESAMPLED_IDS),
    output:
        "04_assemblies/{dataset}/assembly_{assembly_id}.fasta"
    shell:
        """
        echo "Executing...\n"
        canu_trim {input} > {output}
        """

#  Add to troubleshoot necat seg faults
#    benchmark:
#        "benchmarks/{sample}.bwa.benchmark.txt"


rule assemble_with_necat:
    input:
        "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(NECAT_RESAMPLED_IDS),
    output:
        assembly_fasta = "04_assemblies/{dataset}/assembly_{assembly_id}.fasta",
#        assembly_directory = temp(directory("03_necat_assembly_temp/{dataset}/{assembly_id}"))
#    resources:
#        necats=1
    log:
        "logs/{dataset}/necat_{assembly_id}.log"
    conda:
        "../envs/necat.yml"
    threads: 4
    shell:
        """
        ORIGINALPATH=$PWD
        TMPPATH=$(mktemp -d)
        echo "$TMPPATH"
        mkdir -p $TMPPATH
        realpath {input} > $TMPPATH/read_list.txt
        cd $TMPPATH
        necat config config.txt > $ORIGINALPATH/{log} 2>&1
        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
        sed -i "s/GENOME_SIZE=/GENOME_SIZE={config[genome_size]}/" config.txt
        sed -i "s/THREADS=4/THREADS={threads}/" config.txt
        necat bridge config.txt >> $ORIGINALPATH/{log} 2>&1
        cp necat/6-bridge_contigs/polished_contigs.fasta $ORIGINALPATH/{output.assembly_fasta}
        cd $ORIGINALPATH
        rm -r $TMPPATH
        """

# This version doesn't use mktemp    
#        """
#        ORIGINALPATH=$PWD
#        mkdir -p {output.assembly_directory}
#        realpath {input} > {output.assembly_directory}/read_list.txt
#        cd {output.assembly_directory}
#        necat config config.txt > $ORIGINALPATH/{log} 2>&1
#        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
#        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
#        sed -i "s/GENOME_SIZE=/GENOME_SIZE={config[genome_size]}/" config.txt
#        sed -i "s/THREADS=4/THREADS={threads}/" config.txt
#        necat bridge config.txt >> $ORIGINALPATH/{log} 2>&1
#        cp necat/6-bridge_contigs/polished_contigs.fasta $ORIGINALPATH/{output.assembly_fasta}
#        """




READ_NUMS = ["1", "2"]

rule assemble_with_unicycler:
    priority: 1
    input:
        long_reads = "03_subsampled_nanopore_reads/{dataset}/sample_{assembly_id}.fastq",
        short_reads = expand("03_subsampled_illumina_reads/{{dataset}}/sample_{{assembly_id}}_R{read_num}.fastq.gz", read_num=READ_NUMS)
    wildcard_constraints:
        assembly_id="|".join(UNICYCLER_RESAMPLED_IDS),
    output:
        assembly_fasta = "04_assemblies/{dataset}/assembly_{assembly_id}.fasta",
        assembly_gfa = "04_assemblies/{dataset}/assembly_{assembly_id}.gfa",
        assembly_directory = directory("03_unicycler_assembly_temp/{dataset}/{assembly_id}")
    log:
        "logs/{dataset}/unicycler_{assembly_id}.log"
    conda:
        "../envs/unicycler.yml"
    threads: 24
    shell:
        """
        unicycler --threads {threads} -1 {input.short_reads[0]} -2 {input.short_reads[1]} -l {input.long_reads} -o {output.assembly_directory} > {log} 2>&1
        cp {output.assembly_directory}/assembly.fasta {output.assembly_fasta}
        cp {output.assembly_directory}/assembly.gfa {output.assembly_gfa}
        """

rule trycycler_cluster:
    input:
        reads = "02_filtered_nanopore_reads/{dataset}.fastq",
        assemblies = expand("04_assemblies/{{dataset}}/assembly_{assembly_id}.fasta", assembly_id=ALL_RESAMPLED_IDS)
    output:
        output_directory = directory("05_trycycler/{dataset}"),
        done_file = "05_trycycler/{dataset}/done"
    log:
        "logs/{dataset}/tricycler_cluster_{dataset}.log"
    conda:
        "../envs/trycycler.yml"
    threads: 16
    shell:
        "trycycler cluster --threads {threads} --assemblies {input.assemblies} --reads {input.reads} --out_dir 05_trycycler/{wildcards.dataset} > {log} 2>&1 && touch {output.done_file}"

