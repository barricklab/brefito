include: "trim-illumina-reads.smk"

# We need to fix the contig names with this polisher
NORMALIZE_ASSEMBLIES_SCRIPT_PATH = os.path.join(workflow.current_basedir, "..", "scripts", "normalize_assembly.py")

rule polish_with_polypolish:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = expand("01_trimmed_illumina_reads/{{sample}}.R{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        path = temp(directory("06_polypolish/{sample}")),
        polished_fasta = "06_polypolish/{sample}.fasta.polished.temp"
    conda:
        "../envs/polypolish.yml"
    log:
        "logs/polis_polypolish_{sample}.log"
    threads: 16
    shell:
      """
      mkdir -p {output.path}
      cp {input.reference} {output.path}/reference.fasta 
      bwa index {output.path}/reference.fasta 2> {log}
      bwa mem -t 16 -a {output.path}/reference.fasta {input.reads[0]} > {output.path}/alignments_1.sam 2>> {log}
      bwa mem -t 16 -a {output.path}/reference.fasta {input.reads[1]} > {output.path}/alignments_2.sam 2>> {log}
      polypolish {output.path}/reference.fasta 06_polypolish/{wildcards.sample}/alignments_1.sam {output.path}/alignments_2.sam > {output.polished_fasta} 2>> {log}
      """    

rule polish_with_polypolish_postprocess_rename_contigs:
    input:
        reference = "assemblies/{sample}.fasta",
        sample = "06_polypolish/{sample}.fasta.polished.temp"
    output:
        "assemblies/{sample}.fasta.polished"
    conda:
        "../envs/biopython.yml"
    shell:
        "{NORMALIZE_ASSEMBLIES_SCRIPT_PATH} -r {input.reference} -i {input.sample} -o {output} -c"
