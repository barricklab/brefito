include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

import glob
import os.path

BRESEQ_OPTIONS = ""
if 'breseq_options' in config.keys():
    BRESEQ_OPTIONS = config['breseq_options']


import glob 
REFERENCE_FILES = glob.glob("references/*")
#print(REFERENCE_FILES)
REFERENCE_FILES = [ (" -r " + x) for x in REFERENCE_FILES]
REFERENCE_ARGUMENT = ' '.join(REFERENCE_FILES)
#print(REFERENCE_ARGUMENT)

# rule evaluate_breseq_with_nanopore_reads:
#     input:
#         reads = "01_trimmed_nanopore_reads/{sample}.trimmed.fastq.gz"
#     output:
#         intermediate_path = temp(directory("breseq/{sample}_nanopore_reads")),
#         output_path = directory("output/{sample}_nanopore_reads")
#     log: 
#         "logs/breseq_nanopore_reads_{sample}.log"
#     conda:
#         "../envs/breseq.yml"
#     threads: 8
#     shell:
#         """
#         breseq {BRESEQ_OPTIONS} --nanopore --long-read-split-length 500 --long-read-distribute-remainder -j {threads} -o {output.intermediate_path} {input.reads} > {log} 2>&1
#         cp -r {output.intermediate_path}/output {output.output_path}
#        """    

AUTOMATIC_BRESEQ_ARGS = {}

def find_available_read_files(wildcards):
    illumina_path = "illumina_reads/{sample}.R*.fastq.gz"
    illumina_files = glob.glob(illumina_path.format(sample=wildcards.sample))
    nanopore_path = "nanopore_reads/{sample}.fastq.gz"
    nanopore_files = glob.glob(nanopore_path.format(sample=wildcards.sample))

    AUTOMATIC_BRESEQ_ARGS[wildcards.sample] = "";
    if (len(nanopore_files) > 0):
        AUTOMATIC_BRESEQ_ARGS[wildcards.sample] = AUTOMATIC_BRESEQ_ARGS[wildcards.sample] + "-x"

    illumina_input_files = [ os.path.join("01_trimmed_illumina_reads", os.path.basename(s)) for s in illumina_files]
    nanopore_input_files = [ os.path.join("01_trimmed_nanopore_reads", os.path.basename(s)) for s in nanopore_files]

    print(illumina_input_files)
    print(nanopore_input_files)

    return illumina_input_files + nanopore_input_files

def get_breseq_args(sample):
    print(AUTOMATIC_BRESEQ_ARGS[sample])
    return AUTOMATIC_BRESEQ_ARGS[sample]

rule predict_mutations_with_breseq:
    input:
        reads = find_available_read_files
    output:
        breseq = directory("breseq/{sample}"),
        output = directory("output/{sample}")
    log: 
        "logs/evaluate_breseq_{sample}.log"
    conda:
        "../envs/breseq.yml"
    params:
        automatic_breseq_args=lambda wildcards: get_breseq_args(wildcards.sample)
    threads: 8
    shell:
        """
        mkdir -p output
        mkdir -p breseq
        breseq -j {threads} {params.automatic_breseq_args} {BRESEQ_OPTIONS} {REFERENCE_ARGUMENT} -o {output.breseq}  {input.reads} > {log} 2>&1
        cp -r {output.breseq}/output {output.output}
        """    


        # """
        # mkdir -p output
        # mkdir -p breseq
        # BRESEQ_TMP_DIR=$(mktemp -d)
        # echo "$BRESEQ_TMP_DIR"
        # breseq -j {threads} {params.automatic_breseq_args} {BRESEQ_OPTIONS} {REFERENCE_ARGUMENT} -o $BRESEQ_TMP_DIR  {input.reads} > {log} 2>&1
        # # Add the sample name as the TITLE in the GD file - so we can do meaningful gdtools COMPARE
        # #sed -i -e 's/^\(#=GENOME_DIFF.*\)$/\1\n#=TITLE\t{wildcards.sample}/' $BRESEQ_TMP_DIR/output/output.gd
        # #echo "cp -r $BRESEQ_TMP_DIR  {output.breseq}"
        # cp -r $BRESEQ_TMP_DIR {output.breseq}
        # #echo "cp -r {output.breseq}/output {output.output}"
        # cp -r {output.breseq}/output {output.output}
        # """    