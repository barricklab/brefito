try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule merge_all_reads:
    input:
        ["merged_references/{}.fasta".format(s) for s in sample_info.get_sample_list()]
    default_target: True

# This merges all reference files for each sample into one FASTA file and one 
# feature-only GFF3 file (with no nucleotide sequence included).
#
# Is is useful when a tool requires one input files in FASTA format.
# and when you want a GFF feature file for visualization (e.g., for IGV).

rule merge_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = "merged_references/{sample}.fasta",
        gff3 = "merged_references/{sample}.gff3"
    conda:
        "../envs/breseq.yml"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample)
    threads: 1
    shell:
        """
        mkdir -p merged_references
        breseq CONVERT-REFERENCE -o {output.fasta} -f FASTA {params.reference_arguments}
        breseq CONVERT-REFERENCE -o {output.gff3} --no-sequence -f GFF3 {params.reference_arguments} 
        """