try: sample_info
except NameError:
    include: "load-sample-info.smk"

# Does all steps through assembly and clustering

include: "autocycler-subsample-nanopore-reads.smk"

READ_SUBSET_IDS=["01","02","03","04"] # autocycler subsample creates 4 subsampled fastqs for a given fastq file
ASSEMBLERS=["canu","flye","miniasm","necat","nextdenovo","raven"]

rule all_targets_autocycler_cluster:
    input:
        ["autocycler/" + s + "/clustering/" for s in sample_info.get_sample_list()]
    default_target: True

rule canu_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
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
        output_prefix="intermediate-assemblies/{dataset}/canu_assembly_{assembly_id}"
    shell:
        """
        canu.sh {input} {params} {threads} {config[genome_size]}
        """

rule flye_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
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
        output_prefix="intermediate-assemblies/{dataset}/flye_assembly_{assembly_id}"
    shell:
        """
        flye.sh {input} {params} {threads} {config[genome_size]}
        """

rule miniasm_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
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
        output_prefix="intermediate-assemblies/{dataset}/miniasm_assembly_{assembly_id}"
    shell:
        """
        miniasm.sh {input} {params} {threads} {config[genome_size]}
        """

rule necat_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
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
        output_prefix="intermediate-assemblies/{dataset}/necat_assembly_{assembly_id}"
    shell:
        """
        necat.sh {input} {params} {threads} {config[genome_size]}
        """

rule nextdenovo_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
    wildcard_constraints:
        assembly_id="|".join(READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"
    output:
        "intermediate-assemblies/{dataset}/nextdenovo_assembly_{assembly_id}.fasta"
    threads: 16
    log:
        "logs/{dataset}/nextdenovo_{assembly_id}.log"
    params:
        output_prefix="intermediate-assemblies/{dataset}/nextdenovo_assembly_{assembly_id}"
    shell:
        """
        nextdenovo.sh {input} {params} {threads} {config[genome_size]}
        """

rule raven_assemble:
    input:
        "autocycler-nanopore-reads-subsampled/{dataset}/sample_{assembly_id}.fastq"
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
        output_prefix="intermediate-assemblies/{dataset}/raven_assembly_{assembly_id}"
    shell:
        """
        raven.sh {input} {params} {threads} {config[genome_size]}
        """

# compress all of the input assemblies (6 assemblers x 4 subsets each = 24 assemblies per sample) with autocycler

rule autocycler_compress:
    input:
        expand("intermediate-assemblies/{dataset}/{assembler}_assembly_{assembly_id}.fasta",
               dataset=sample_info.get_sample_list(),
               assembly_id=READ_SUBSET_IDS,
               assembler=ASSEMBLERS)
    output:
        "autocycler/{dataset}/input_assemblies.gfa"
    params:
        input_assemblies="intermediate-assemblies/{dataset}/",
        output_directory="autocycler/{dataset}/"
    log:
        "logs/{dataset}/autocycler_compress_{dataset}.log"
    threads: 16
    shell:
        """
        autocycler compress -i {params.input_assemblies} -a {params.output_directory}
        """

rule autocycler_cluster:
    input:
        "autocycler/{dataset}/input_assemblies.gfa"
    output:
        "autocycler/{dataset}/clustering/"
    params:
        input_directory="autocycler/{dataset}/"
    log:
        "logs/{dataset}/autocycler_cluster_{dataset}.log"
    threads: 16
    shell:
        """
        autocycler cluster -a {params.input_directory}
        """
