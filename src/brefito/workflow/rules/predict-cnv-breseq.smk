try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "predict-mutations-breseq.smk"

import os.path

rule all_predict_cnv_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/cnv/" + s + "/CNV_csv/coverage_break_pts.csv" for s in sample_info.get_sample_list()]
    default_target: True

rule predict_cnv_breseq:
    input:
        "breseq-" + sample_info.get_reference_prefix() + "/cnv/{sample}/coverage.tsv",
    output:
        "breseq-" + sample_info.get_reference_prefix() + "/cnv/{sample}/CNV_csv/coverage_break_pts.csv"
    params:
        path =  "breseq-" + sample_info.get_reference_prefix() + "/cnv/{sample}"
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}-predict-cnv.log"
    conda:
        "../envs/breseq-ext-cnv.yml"
    threads: 1
    shell:
        """
        breseq-ext-cnv -i {input} -o {params.path} -w 100 -s 100 > {log} 2>&1
        """

rule create_coverage_file_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
    output:
        temp("breseq-" + sample_info.get_reference_prefix() + "/cnv/{sample}/coverage.tsv")
    params:
        path =  "breseq-" + sample_info.get_reference_prefix() + "/cnv/{sample}"
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}-create-coverage-file.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq BAM2COV --table --resolution 0 --bam {input.bam} --fasta {input.fasta} --output {params.path} > {log} 2>&1
        mv {params.path}/*.tab {output}
        """
