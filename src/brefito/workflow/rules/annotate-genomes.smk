# There are some bugs generating Genbank files
# with versions of breseq < 0.38.3, 
# which is why this makes GFF3 only right now!

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_annotate_genomes:
    input:
        ["annotated_" + sample_info.get_reference_prefix() + "/" + s + ".gff3" for s in sample_info.get_sample_list()]
    default_target: True

#Convert the references to fasta format (so we can re-annotate if needed)
rule convert_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = temp("annotated_" + sample_info.get_reference_prefix() + "_fasta/{sample}.fasta")
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_convert_references.log"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample),
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -o {output.fasta} -f FASTA {params.reference_arguments} > {log} 2>&1
        """

rule annotate_with_prokka:
    input:
        "annotated_" + sample_info.get_reference_prefix() + "_fasta/{sample}.fasta"
    output:
        dir = temp(directory("annotated_" + sample_info.get_reference_prefix() + "_prokka/{sample}")),
        prokka_annotated_reference = temp("annotated_" + sample_info.get_reference_prefix() + "_prokka/{sample}/reference.gff")
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
        "annotated_" + sample_info.get_reference_prefix() + "_fasta/{sample}.fasta"
    output:
        dir = temp(directory("annotated_" + sample_info.get_reference_prefix() + "_isescan/{sample}")),
        isescan_annotated_csv = temp("annotated_" + sample_info.get_reference_prefix() + "_isescan/{sample}/" +  "annotated_" + sample_info.get_reference_prefix() + "_fasta" + "/{sample}.fasta.csv")
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_isescan.log"
    conda:
        "../envs/isescan.yml"
    threads: 8
    shell:
        """
        isescan.py --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1
        """

#       isescan.py --removeShortIS --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1


rule combine_annotation_with_breseq:
    input:
        prokka = "annotated_" + sample_info.get_reference_prefix() + "_prokka/{sample}/reference.gff",
        isescan = "annotated_" + sample_info.get_reference_prefix() + "_isescan/{sample}/" +  "annotated_" + sample_info.get_reference_prefix() + "_fasta" + "/{sample}.fasta.csv",
        prokka_dir = "annotated_" + sample_info.get_reference_prefix() + "_prokka/{sample}",
        isescan_dir = "annotated_" + sample_info.get_reference_prefix() + "_isescan/{sample}"
    output:
        "annotated_" + sample_info.get_reference_prefix() + "/{sample}.gff3"
    log:
        "logs/" + "annotated_" + sample_info.get_reference_prefix() + "_{sample}_combine_annotation_with_breseq.log"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f GFF3 -s {input.isescan} -o {output} {input.prokka} > {log} 2>&1
        """
