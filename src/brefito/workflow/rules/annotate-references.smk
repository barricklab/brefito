rule annotate_with_prokka:
    input:
        "references/{sample}.fasta"
    output:
        dir = directory("intermediate/prokka/{sample}"),
        prokka_annotated_reference = "intermediate/prokka/{sample}/reference.gbk"
    conda:
        "../envs/prokka.yml"
    threads: 32
    shell:
        """
        prokka --cpus {threads} --prefix reference --force --outdir {output.dir} {input}
        """

rule annotate_with_isescan:
    input:
        "references/{sample}.fasta"
    output:
        dir = directory("intermediate/isescan/{sample}"),
        isescan_annotated_csv = "intermediate/isescan/{sample}/references/{sample}.fasta.csv"
    conda:
        "../envs/isescan.yml"
    threads: 32
    shell:
        """
        isescan.py --removeShortIS --nthread {threads} --seqfile {input} --output {output.dir}
        """

rule combine_annotation_with_breseq:
    input:
        prokka = "intermediate/prokka/{sample}/reference.gbk",
        isescan = "intermediate/isescan/{sample}/references/{sample}.fasta.csv"
    output:
        "output/annotated_references/{sample}.gbk"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f GENBANK -s {input.isescan} -o {output} {input.prokka}  
        """
