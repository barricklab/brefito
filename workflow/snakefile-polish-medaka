import glob
import os.path


input_files=glob.glob("input/*.long_reads.fastq.gz")
#print(input_files)
DATASETS=[]
print("Running on the following datasets...")
for this_input_file in input_files:
    this_base_name=os.path.basename(this_input_file)
    #print(this_base_name)
    res = re.findall("(.+).long_reads.fastq.gz", this_base_name)
    if not res: continue
    print("  ", res[0])
    DATASETS.append(res[0])


rule all:
  input:
        expand("output/{dataset}.fasta", dataset=DATASETS)

rule polish_with_medaka:
    input:
        reference = "05_trycycler/{sample}/cluster_{cluster}/7_final_consensus.fasta",
        reads = "01_trimmed_reads/{sample}_long_reads.fastq"
    output:
        "05_trycycler/{sample}/cluster_{cluster}/8_medaka.fasta"
    conda:
        "envs/medaka.yml"
    threads: 16
    shell:
      """
      medaka_consensus -i {input.reads} -d {input.reference} -o 05_trycycler/{wildcards.sample}/cluster_{wildcards.cluster}/medaka -m r941_min_sup_g507 -t {threads}
      mv 05_trycycler/{wildcards.sample}/cluster_{wildcards.cluster}/medaka/consensus.fasta {output}
      rm -r 05_trycycler/{wildcards.sample}/cluster_{wildcards.cluster}/medaka 
      rm 05_trycycler/{wildcards.sample}/cluster_{wildcards.cluster}/*.fai
      rm 05_trycycler/{wildcards.sample}/cluster_{wildcards.cluster}/*.mmi
      """    
   
def find_clusters(wildcards):
    clusters = [str(x) for x in glob.glob(f"05_trycycler/{wildcards.dataset}/cluster_*/")]
    clusters = list(map(lambda x: x + "8_medaka.fasta", clusters))
    print(clusters)
    return sorted(clusters) 

rule combine_medaka_polished_contigs:
    input:
        polished_contigs = find_clusters
    output:
        "output/{dataset}.fasta"
    shell:
        "cat {input.polished_contigs} > {output}"

