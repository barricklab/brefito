include: "trim-illumina-reads.smk"

@@@ BROKEN - needs to merge trimmed illumina reads
@@@ Should also be converted to use pypolca!

rule polish_with_polca:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = expand("illumina-reads-trimmed/{{sample}}.R{read_num}.fastq.gz", read_num=["1", "2"])
    output:
        "assemblies/{sample}.fasta.polished"
    log:
        "logs/polish_polca_{sample}.log"
    conda:
        "../envs/masurca.yml"
    threads: 16
    shell:
      """
      MAIN_DIR=$PWD
      WORK_DIR=`mktemp -d`
      cd $WORK_DIR
      polca.sh -a $MAIN_DIR/{input.reference} -r "$MAIN_DIR/{input.reads[0]} $MAIN_DIR/{input.reads[1]} " -t 16 -m 1G
      rm  $MAIN_DIR/{input.reference}.fai
      mv *.PolcaCorrected.fa $MAIN_DIR/{output}
      """    
