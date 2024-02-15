try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_mutate_genomes_gdtools:
    input:
        ["mutants/" + s + ".gff3" for s in sample_info.get_sample_list()]
    default_target: True

rule mutate_genomes_gdtools:
    input:
        genomediff = "genome-diffs/{sample}.gd",
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "mutants/{sample}.gff3"
    log: 
        "logs/mutate-genomes-gdtools-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: 8
    shell:
        """
        gdtools APPLY -o {output} {params.reference_arguments} -f GFF3 {input.genomediff} -r {input.references} > {log} 2>&1
        """
