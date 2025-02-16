try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "download-data.smk"

#These control how the contigs are named after normalization
# --config RENAME_CONTIGS=sample/reference

print("=========================================================")
print("Config options for workflow: normalize-assembles")
print("=========================================================")
print("--config RENAME_CONTIGS=<options>")
print("         SAMPLE renames contigs to the name of the sample,")
print("           of sample_# if there is more than one contig.")
print("         REFERENCE renames contigs to the name of the")
print("           corresponding contig in the reference file")
print("           (assuming they are the same sorted by length")
print("         The default is to not change contig names.")
print("=========================================================")

NORMALIZE_NAMING_OPTIONS = ""
RENAME_CONTIGS = ""
if 'RENAME_CONTIGS' in brefito_config.keys():
    RENAME_CONTIGS = brefito_config['RENAME_CONTIGS'].upper()

ruleorder: normalize_genomes_no_rename > normalize_genomes_rename_sample > normalize_genomes_rename_reference
if RENAME_CONTIGS == "SAMPLE":
    print("Renaming contigs by sample.")
    ruleorder: normalize_genomes_rename_sample > normalize_genomes_no_rename > normalize_genomes_rename_reference
elif RENAME_CONTIGS == "REFERENCE":
    print("Renaming contigs by reference.")
    ruleorder: normalize_genomes_rename_reference > normalize_genomes_rename_sample > normalize_genomes_no_rename

rule all_normalize_genomes:
    input:
        ["assemblies/" + s + ".fasta.normalized" for s in sample_info.get_sample_list() ]
    default_target: True
    params:
        reference_path = "normalize-merged-references"
    shell:
        "rm -rf {params.reference_path}"

rule create_merged_reference_fasta:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "normalize-merged-references/{sample}.fasta"
    threads: 1
    log: 
        "logs/normalize-merged-references-{sample}.log"
    params:
        reference_path = "normalize-merged-references"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        mkdir -p {params.reference_path}
        breseq CONVERT-REFERENCE -f FASTA -o {output} {input.references} > {log} 2>&1
        """

rule normalize_genomes_rename_sample:
    input:
        reference = "normalize-merged-references/{sample}.fasta",
        assembly = "assemblies/{sample}.fasta"
    output:
        "assemblies/{sample}.fasta.normalized"
    threads: 1
    log: 
        "logs/normalize-assemblies-{sample}.log"
    shell:
        """
        normalize_assembly -r {input.reference} -i {input.assembly} -o {output} -s -x -n {wildcards.sample} > {log} 2>&1
        """

rule normalize_genomes_rename_reference:
    input:
        reference = "normalize-merged-references/{sample}.fasta",
        assembly = "assemblies/{sample}.fasta"
    output:
        "assemblies/{sample}.fasta.normalized"
    threads: 1
    log: 
        "logs/normalize-assemblies-{sample}.log"
    shell:
        """
        normalize_assembly -r {input.reference} -i {input.assembly} -o {output} -s -x -c > {log} 2>&1
        """

rule normalize_genomes_no_rename:
    input:
        reference = "normalize-merged-references/{sample}.fasta",
        assembly = "assemblies/{sample}.fasta"
    output:
        "assemblies/{sample}.fasta.normalized"
    threads: 1
    log: 
        "logs/normalize-assemblies-{sample}.log"
    shell:
        """
        normalize_assembly -r {input.reference} -i {input.assembly} -o {output} -s -x > {log} 2>&1
        """