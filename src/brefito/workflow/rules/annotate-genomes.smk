try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_annotate_genomes:
    input:
        ["annotated_" + sample_info.get_reference_prefix() + s + ".gbk" for s in sample_info.get_sample_list()]
    default_target: True

#Convert the references to fasta format (so we can re-annotate if needed)
rule convert_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = temp("annotated_" + sample_info.get_reference_prefix() + "/{sample}.fasta")
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample)
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -o {output.fasta} -f FASTA {params.reference_arguments}
        """

rule annotate_with_prokka:
    input:
        "references/{sample}.fasta"
    output:
        dir = temp(directory("annotated_" + sample_info.get_reference_prefix() + "/prokka/{sample}")),
        prokka_annotated_reference = temp("annotated_" + sample_info.get_reference_prefix() + "/prokka/{sample}/reference.gbk")
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_prokka.log"
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
        dir = temp(directory("annotated_" + sample_info.get_reference_prefix() + "/isescan/{sample}")),
        isescan_annotated_csv = temp("annotated_" + sample_info.get_reference_prefix() + "/isescan/{sample}/references/{sample}.fasta.csv")
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_isescan.log"
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
        prokka = "annotated_" + sample_info.get_reference_prefix() + "/prokka/{sample}/reference.gbk",
        isescan = "annotated_" + sample_info.get_reference_prefix() + "isescan/{sample}/references/{sample}.fasta.csv"
    output:
        "annotated_" + sample_info.get_reference_prefix() + "/{sample}.gbk"
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_combine_annotation_with_breseq.log"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f GENBANK -s {input.isescan} -o {output} {input.prokka} > {log} 2>&1
        """
