rule annotate_with_prokka:
    input:
        "references/{sample}.fasta"
    output:
        dir = directory("intermediate/prokka/{sample}"),
        prokka_annotated_reference = "intermediate/prokka/{sample}/reference.gbk"
    log:
        "logs/{sample}/prokka.log"
    conda:
        "../envs/prokka.yml"
    threads: 4
    shell:
        """
        prokka --cpus {threads} --prefix reference --force --outdir {output.dir} {input} > {log} 2>&1
        """

rule annotate_with_isescan:
    input:
        "references/{sample}.fasta"
    output:
        dir = directory("intermediate/isescan/{sample}"),
        isescan_annotated_csv = "intermediate/isescan/{sample}/references/{sample}.fasta.csv"
    log:
        "logs/{sample}/isescan.log"
    conda:
        "../envs/isescan.yml"
    threads: 4
    shell:
        """
        isescan.py --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1
        """

#       isescan.py --removeShortIS --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1


rule combine_annotation_with_breseq:
    input:
        prokka = "intermediate/prokka/{sample}/reference.gbk",
        isescan = "intermediate/isescan/{sample}/references/{sample}.fasta.csv"
    output:
        "output/annotated_references/{sample}.gbk"
    log:
        "logs/{sample}/combine_annotation_with_breseq.log"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f GENBANK -s {input.isescan} -o {output} {input.prokka} > {log} 2>&1
        """
