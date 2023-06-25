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

rule polish_with_polca:
    input:
        reference = "input/{dataset}.fasta",
        reads = expand("01_trimmed_reads/{{dataset}}.short_reads_{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        "output/{dataset}.fasta"
    conda:
        "envs/masurca.yml"
    threads: 16
    shell:
      """
      MAIN_DIR=$PWD
      WORK_DIR=`mktemp -d`
      cd $WORK_DIR
      polca.sh -a $MAIN_DIR/{input.reference} -r "$MAIN_DIR/{input.reads[0]} $MAIN_DIR/{input.reads[1]} " -t 16 -m 1G
      mv *.PolcaCorrected.fa $MAIN_DIR/{output}
      """    
