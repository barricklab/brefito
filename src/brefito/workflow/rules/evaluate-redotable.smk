if not "redotable_path" in config.keys():
    print("You must supply a path to redotoable!\n Like this: --config redotable_path=/path/to/redotable\n")

REDOTABLE_PATH = config["redotable_path"]

rule evaluate_redotable:
    input:
        assembly = "assemblies/{sample}.fasta",
        reference = "reference.fasta"
    output:
        plot = "evaluate/dot_plot/{sample}.svg"
    log:
        "logs/{sample}/redotable.log"
#    No conda package, so we pass path as config variable
#    conda:
#        "../envs/redotable.yml"
    threads: 1
    shell:
        "perl {REDOTABLE_PATH}/redotable --width 2000 --height 2000 --svg {input.reference} {input.assembly} {output.plot} > {log} 2>&1"