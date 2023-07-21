rule normalize_assemblies:
    input:
        reference = "reference.fasta",
        assembly = "assemblies/{sample}.fasta"
    output:
        "assemblies/{sample}.fasta.normalized"
    shell:
        "normalize_assembly -r {input.reference} -i {input.assembly} -o {output} -s -c -x"
