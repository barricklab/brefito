try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_validate_genomediff_gdtools:
    input:
        ["genome-diffs-validate/" + s + ".txt" for s in sample_info.get_sample_list()]
    default_target: True

rule _validate_genomediff_gdtools:
    input:
        genomediff = "genome-diffs/{sample}.gd",
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "genome-diffs-validate/{sample}.txt"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        gdtools VALIDATE -r {input.references} {input.genomediff} > {output} 2>&1 
        """
