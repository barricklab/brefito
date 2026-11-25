try: sample_info
except NameError:
    include: "load-sample-info.smk"

# Does all steps for autocycler after assembly

include: "autocycler-subsample-nanopore-reads.smk"

READ_SUBSET_IDS=["01","02","03","04"] # autocycler subsample creates 4 subsampled fastqs for a given fastq file
ASSEMBLERS=["canu","flye","miniasm","necat","metamdbg","raven"]

# Assemblers run for a dataset, by its long-read type: NECAT assembles nanopore reads only, and
# metaMDBG only PacBio HiFi and high-accuracy (nanopore-hq) nanopore reads.
def get_assemblers(dataset):
    long_read_type = sample_info.get_long_read_type(dataset)
    return [a for a in ASSEMBLERS
            if (a != "necat" or long_read_setting(long_read_type, "necat"))
            and (a != "metamdbg" or long_read_setting(long_read_type, "metamdbg"))]

# autocycler helper --read_type for the dataset's reads (ont_r9 for nanopore and nanopore-raw,
# ont_r10 for nanopore-hq, pacbio_clr, pacbio_hifi)
def get_autocycler_read_type_argument(dataset):
    return "--read_type " + long_read_setting(sample_info.get_long_read_type(dataset), "autocycler")

MAX_CONTIGS_ARGUMENT="--max_contigs 40"

rule all_targets_autocycler:
    input:
        ["assemblies/" + s + ".fasta" for s in sample_info.get_samples_with_long_reads()]
    default_target: True

rule canu_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/canu_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/canu_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/canu_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper canu --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

rule flye_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/flye_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/flye_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/flye_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper flye --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

rule miniasm_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/miniasm_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/miniasm_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/miniasm_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper miniasm --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

rule necat_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/necat_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/necat_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/necat_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper necat --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

#rule nextdenovo_assemble:
#    input:
#        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
#    wildcard_constraints:
#        assembly_id="|".join(READ_SUBSET_IDS)
#    conda:
#        "../envs/autocycler.yml"
#    output:
#        "intermediate-assemblies/{dataset}/nextdenovo_assembly_{assembly_id}.fasta"
#    threads: 16
#    log:
#        "logs/{dataset}/nextdenovo_{assembly_id}.log"
#    params:
#        output_prefix="intermediate-assemblies/{dataset}/nextdenovo_assembly_{assembly_id}"
#    shell:
#        """
#        autocycler helper nextdenovo --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {config[genome_size]} > {log} 2>&1
#        """

rule metamdbg_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/metamdbg_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/metamdbg_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/metamdbg_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper metamdbg --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

rule raven_assemble:
    input:
        lambda wildcards: get_autocycler_subsampled_reads(wildcards.dataset, wildcards.assembly_id)
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/raven_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/raven_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/raven_assembly_{assembly_id}",
        genome_size=lambda wildcards: get_genome_size(wildcards.dataset),
        read_type=lambda wildcards: get_autocycler_read_type_argument(wildcards.dataset)
    shell:
        """
        autocycler helper raven --reads {input} --out_prefix {params.output_prefix} --threads {threads} --genome_size {params.genome_size} {params.read_type} > {log} 2>&1
        """

# compress all of the input assemblies (6 assemblers x 4 subsets each = 24 assemblies per sample) with autocycler

rule autocycler_compress:
    input:
        # only this dataset's assemblies, from the assemblers its read type allows
        lambda wildcards: expand("intermediate-assemblies/{dataset}/{assembler}_assembly_{assembly_id}.fasta",
               dataset=wildcards.dataset,
               assembly_id=READ_SUBSET_IDS,
               assembler=get_assemblers(wildcards.dataset))
    output:
        "autocycler/{dataset}/input_assemblies.gfa"
    conda:
        "../envs/autocycler.yml"
    params:
        input_assemblies="intermediate-assemblies/{dataset}/",
        output_directory="autocycler/{dataset}/"
    log:
        "logs/{dataset}/autocycler_compress_{dataset}.log"
    threads: 16
    shell:
        """
        autocycler compress --threads {threads} {MAX_CONTIGS_ARGUMENT} -i {params.input_assemblies} -a {params.output_directory} > {log} 2>&1
        """

rule autocycler_all_steps:
    input:
        input_gfa = "autocycler/{dataset}/input_assemblies.gfa"
    output:
        "assemblies/{dataset}.fasta"
    conda:
        "../envs/autocycler.yml"
    log:
        "logs/{dataset}/autocycler.log"
    shell:
       """
       sample={wildcards.dataset}
       input_directory="autocycler/$sample"

       # run cluster
       autocycler cluster {MAX_CONTIGS_ARGUMENT} -a $input_directory > {log} 2>&1

       for c in "$input_directory/clustering/qc_pass/cluster_"*; do
           autocycler trim  --threads {threads} -c "$c" > {log} 2>&1
           #autocycler dotplot -i "$c"/1_untrimmed.gfa -o "$c"/1_untrimmed.png > {log} 2>&1
           autocycler dotplot -i "$c"/2_trimmed.gfa -o "$c"/2_trimmed.png > {log} 2>&1
           autocycler resolve -c "$c" > {log} 2>&1
       done

       output_directory="autocycler/$sample/output"

       if [ ! -d "$output_directory" ]; then
           mkdir -p "$output_directory"
       fi

       autocycler combine -a "$output_directory" -i "$input_directory"/clustering/qc_pass/cluster_*/5_final.gfa > {log} 2>&1 

       mkdir -p assemblies
       cp $output_directory/consensus_assembly.fasta assemblies/{wildcards.dataset}.fasta
       cp $output_directory/consensus_assembly.gfa assemblies/{wildcards.dataset}.gfa
       cp $output_directory/consensus_assembly.yaml assemblies/{wildcards.dataset}.yaml
       """
