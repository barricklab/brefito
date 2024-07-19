try: sample_info
except NameError: 
    include: "load-sample-info.smk"

from os import system
system("rm blast-" + sample_info.get_reference_prefix() + "/blast_output*")

"blast-" + sample_info.get_reference_prefix() + "/blast_output.txt"

INPUT = brefito_config['INPUT'] 

rule create_merged_reference_fasta:
    input:
        sample_info.get_all_reference_list()
    output:
        "blast-" + sample_info.get_reference_prefix() + "/merged_references.fasta"
    threads: 1
    log: 
        "logs/merge-references-" + sample_info.get_reference_prefix() + ".log"
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -f FASTA -o {output} {input}
        """

rule references_format_blastdb:
    input:
        "blast-" + sample_info.get_reference_prefix() + "/merged_references.fasta"
    output:
        "blast-" + sample_info.get_reference_prefix() + "/blastdb.ndb"
    threads: 1
    log: 
        "logs/blast-formatdb-" + sample_info.get_reference_prefix() + ".log"
    conda:
        "../envs/blast.yml"
    params:
        blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb"
    shell:
        """
        makeblastdb -dbtype nucl -in {input} -out {params.blastdb}
        """

rule search_references_blast:
    default_target: True
    input:
       blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb.ndb"
    output:
        format1 = "blast-" + sample_info.get_reference_prefix() + "/blast_output.txt",
        format7 = "blast-" + sample_info.get_reference_prefix() + "/blast_output.tsv"
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + ".log"
    conda:
        "../envs/blast.yml"
    threads: 12
    params:
        blastdb = "blast-" + sample_info.get_reference_prefix() + "/blastdb",
        input = "blast-" + sample_info.get_reference_prefix() + "/blast_input.txt"
    shell:
        """
        echo "{INPUT}" > {params.input}
        blastn -db {params.blastdb} -query {params.input} -out {output.format1}
        cat {output.format1}
        blastn -db {params.blastdb} -outfmt 7 -query {params.input} -out {output.format7}
        cat {output.format7}
        """