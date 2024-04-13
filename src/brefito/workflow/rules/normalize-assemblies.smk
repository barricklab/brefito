try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_normalize_assemblies:
    input:
        ["assemblies/" + s + ".fasta.normalized" for s in sample_info.get_sample_list() ]
    default_target: True

rule normalize_assemblies:
    input:
        reference = "reference.fasta",
        assembly = "assemblies/{sample}.fasta"
    output:
        "assemblies/{sample}.fasta.normalized"
    shell:
        "normalize_assembly -r {input.reference} -i {input.assembly} -o {output} -s -c -x -n {wildcards.sample}"
