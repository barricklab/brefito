rule evaluate_inspector:
    input:
        contigs = "assemblies/{sample}.fasta",
        reads = "01_trimmed_nanopore_reads/{sample}.fastq",
        reference = "reference.fasta"
    output:
        dir = directory("inspector_assembly_evaluation/{sample}"),
        file = "inspector_assembly_evaluation/{sample}/summary_statistics" 
    log:
        "logs/{sample}/inspector.log"
    conda:
        "../envs/inspector.yml"
    threads: 12
    shell:
        "inspector.py -t 24 --skip_base_error -c {input.contigs} -r {input.reads} -o {output.dir} --ref {input.reference} > {log} 2>&1"