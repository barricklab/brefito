try: sample_info
except NameError: 
    include: "load-sample-info.smk"

#from os import system
#system("rm blast-" + sample_info.get_reference_prefix() + "/blast_output*")

QUERY_FILE_PATH = brefito_config['QUERY_FILE_PATH'] 

rule create_merged_reference_fasta:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_sequences.fasta"
    threads: 1
    log: 
        "logs/merge-references-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        blastdb_path = "blast-" + sample_info.get_reference_prefix() + "/blastdb"
    shell:
        """
        mkdir -p {params.blastdb_path}
        breseq CONVERT-REFERENCE -f FASTA -o {output} {input}
        """

rule references_format_blastdb:
    input:
        "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_sequences.fasta"
    output:
        "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_blastdb.ndb"
    threads: 1
    log: 
        "logs/blast-formatdb-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/blast.yml"
    params:
        blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_blastdb"
    shell:
        """
        makeblastdb -dbtype nucl -in {input} -out {params.blastdb}
        """

rule search_references_blast:
    input:
        query = QUERY_FILE_PATH,
        blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_blastdb.ndb"
    output:
        format1 = "blast-" + sample_info.get_reference_prefix() + "/blast_{sample}_output.html",
        format7 = "blast-" + sample_info.get_reference_prefix() + "/blast_{sample}_output.tsv"
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/blast.yml"
    threads: 12
    params:
        blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb/{sample}_blastdb",
    shell:
        """
        blastn -db {params.blastdb} -html -query {input.query} -out {output.format1}
        blastn -db {params.blastdb} -outfmt 7 -query {input.query} -out {output.format7}
        cat {output.format7}
        """

rule all_search_references_blast:
    input:
        ["blast-" + sample_info.get_reference_prefix() + "/blast_" + s + "_output.html" for s in sample_info.get_sample_list()]
    default_target: True