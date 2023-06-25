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

rule map_short_reads_for_pilon:
    input:
        reference = "input/{dataset}.fasta",
        short_reads = expand("01_trimmed_reads/{{dataset}}.short_reads_{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        "06_pilon/{dataset}/mapped_short_reads.sam"
    conda:
        "envs/pilon.yml"
    threads: 16
    shell:
        "bowtie2-build {input.reference} 06_pilon/{wildcards.dataset}/reference && bowtie2 -X 1000 -p {threads} -x 06_pilon/{wildcards.dataset}/reference -1 {input.short_reads[0]} -2 {input.short_reads[1]} -S {output}"

rule map_long_reads_for_pilon:
    input:
       reference = "input/{dataset}.fasta",
       long_reads = "01_trimmed_reads/{dataset}.long_reads.fastq"
    output:
       "06_pilon/{dataset}/mapped_long_reads.sam"
    conda:
        "envs/pilon.yml"
    threads: 16
    shell:
       "minimap2 -ax map-ont -t {threads} {input.reference} {input.long_reads} > {output}"

rule sam_to_bam_for_long_reads_for_pilon:
     input:
        "06_pilon/{dataset}/{mapped}.sam"
     output:
        "06_pilon/{dataset}/{mapped}.bam"     
     conda:
        "envs/pilon.yml" 
     threads: 16  
     shell:
        "samtools view -bS --threads {threads} {input} > {output}"

rule sort_bam_for_pilon:
     input:
        "06_pilon/{dataset}/{mapped}.bam"
     output:
        "06_pilon/{dataset}/{mapped}.sorted.bam"
     conda:
        "envs/pilon.yml"
     threads: 16
     shell:
        "samtools sort -@ {threads} -o {output} {input}"  

rule index_bam_for_pilon:
     input:
        "06_pilon/{dataset}/{mapped}.sorted.bam"
     output:
        "06_pilon/{dataset}/{mapped}.sorted.bam.bai"
     conda:
        "envs/pilon.yml"
     threads: 16
     shell:
        "samtools index -@ {threads} {input}" 

rule polish_with_pilon:
    input:
        reference = "input/{dataset}.fasta",
        short_reads_bam = "06_pilon/{dataset}/mapped_short_reads.sorted.bam",
        long_reads_bam = "06_pilon/{dataset}/mapped_long_reads.sorted.bam",
        short_reads_bam_index = "06_pilon/{dataset}/mapped_short_reads.sorted.bam.bai",
        long_reads_bam_index = "06_pilon/{dataset}/mapped_long_reads.sorted.bam.bai"
    output:
        "output/{dataset}.fasta"
    conda:
        "envs/pilon.yml"
    threads: 16
    shell:
      "pilon -Xmx100G --threads {threads} --outdir 06_pilon/{wildcards.dataset} --genome {input.reference} --frags {input.short_reads_bam} --fix all"
