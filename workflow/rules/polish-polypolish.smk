include: "snakefile-trim-short-reads"

import glob
import os.path


input_files=glob.glob("input/*.fasta")
#print(input_files)
DATASETS=[]
print("Running on the following datasets...")
for this_input_file in input_files:
    this_base_name=os.path.basename(this_input_file)
    #print(this_base_name)
    res = re.findall("(.+).fasta", this_base_name)
    if not res: continue
    print("  ", res[0])
    DATASETS.append(res[0])



rule all:
  input:
        expand("output/{dataset}.fasta", dataset=DATASETS)

rule polish_with_polypolish:
    input:
        reference = "input/{dataset}.fasta",
        reads = expand("01_trimmed_reads/{{dataset}}.short_reads_{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        "output/{dataset}.fasta"
    conda:
        "envs/polypolish.yml"
    threads: 16
    shell:
      """
      mkdir -p 06_polypolish/{wildcards.dataset}
      mkdir -p output
      cp {input.reference} 06_polypolish/{wildcards.dataset}/reference.fasta 
      bwa index 06_polypolish/{wildcards.dataset}/reference.fasta
      bwa mem -t 16 -a 06_polypolish/{wildcards.dataset}/reference.fasta {input.reads[0]} > 06_polypolish/{wildcards.dataset}/alignments_1.sam
      bwa mem -t 16 -a 06_polypolish/{wildcards.dataset}/reference.fasta {input.reads[1]} > 06_polypolish/{wildcards.dataset}/alignments_2.sam
      polypolish 06_polypolish/{wildcards.dataset}/reference.fasta 06_polypolish/{wildcards.dataset}/alignments_1.sam 06_polypolish/{wildcards.dataset}/alignments_2.sam > {output}
      """    
