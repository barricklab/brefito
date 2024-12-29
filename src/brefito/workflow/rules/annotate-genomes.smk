# There are some bugs generating Genbank files
# with versions of breseq < 0.38.3, 
# which is why this makes GFF3 only right now!

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule all_annotate_genomes:
    input:
        ["annotated-" + sample_info.get_reference_prefix() + "/" + s + ".gbk" for s in sample_info.get_sample_list()]
    default_target: True

#Convert the references to fasta format (so we can re-annotate if needed)
rule convert_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = temp("annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "_{sample}_convert_references.log"
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
        "annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta"
    output:
        dir = temp(directory("annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}")),
        prokka_annotated_reference = temp("annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}/reference.gff")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-prokka.log"
    conda:
        "../envs/prokka.yml"
    threads: 4
    shell:
        """
        prokka --cpus {threads} --prefix reference --force --outdir {output.dir} {input} > {log} 2>&1
        """

rule annotate_with_isescan:
    input:
        "annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta"
    output:
        dir = temp(directory("annotated-" + sample_info.get_reference_prefix() + "-isescan/{sample}")),
        isescan_annotated_csv = temp("annotated-" + sample_info.get_reference_prefix() + "-isescan/{sample}/" +  "annotated-" + sample_info.get_reference_prefix() + "-fasta" + "/{sample}.fasta.csv")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-isescan.log"
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
        prokka = "annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}/reference.gff",
        isescan = "annotated-" + sample_info.get_reference_prefix() + "-isescan/{sample}/" +  "annotated-" + sample_info.get_reference_prefix() + "-fasta" + "/{sample}.fasta.csv",
        prokka_dir = "annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}",
        isescan_dir = "annotated-" + sample_info.get_reference_prefix() + "-isescan/{sample}"
    output:
        "annotated-" + sample_info.get_reference_prefix() + "/{sample}.gbk"
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-combine-annotation-with-breseq.log"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f GENBANK -s {input.isescan} -o {output} {input.prokka} > {log} 2>&1
        """
