try: sample_info
except NameError: 
    include: "load-sample-info.smk"

#include: "annotate-genomes.smk"

# Default; BRESEQ_OPTIONS can be set per sample in data.csv or with --config
DEFAULT_BRESEQ_OPTIONS = "--max-displayed-reads 40 --junction-score-cutoff 0 --junction-minimum-pos-hash-score 10 --brief-html-output"


rule convert_annotated_references_to_genbank:
    input:
        references = "annotated-" + sample_info.get_reference_prefix() + "/{sample}.gff3"
    output:
        genbank = "annotated-" + sample_info.get_reference_prefix() + "-genbank/{sample}.gbk"
    log: 
        "logs/breseq-convert-references-to-genbank-{sample}.log"
    conda:
        BRESEQ_ENV
    threads: 1
    shell:
        """
        breseq convert-reference -f GENBANK -o {output.genbank} {input.references} > {log} 2>&1
        """

rule extract_mobile_elements:
    input:
        "annotated-" + sample_info.get_reference_prefix() + "-genbank/{sample}.gbk"
    output:
        "mobile-element-genbanks/{sample}.gbk"
    log: 
        "logs/extract-mobile-elements-{sample}.log"
    threads: 1
    shell:
        """
        extract_mobile_elements -t genbank --input {input} --output {output}  > {log} 2>&1
        """


rule all_compare_assemblies_breseq:
    input:
        ["breseq-compare-" + sample_info.get_reference_prefix() + "/html/" + s + "/output.done" for s in sample_info.get_sample_list()]
    default_target: True

# Annotate assemblies

rule simulate_reads_breseq:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        reads_gz = "simulated-reads-" + sample_info.get_reference_prefix() + "/{sample}.fastq.gz"
    log:
        "logs/breseq-simulate-reads-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        reads = "simulated-reads-" + sample_info.get_reference_prefix() + "/{sample}.fastq",
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: 8
    shell:
        """
        breseq simulate-reads -m tiled -l 200 -c 80 {params.reference_arguments} -o {params.reads}  > {log} 2>&1
        gzip {params.reads}
        """


rule predict_mutations_breseq:
    input:
        reads = "simulated-reads-" + sample_info.get_reference_prefix() + "/{sample}.fastq.gz",
        main_reference = sample_info.get_main_reference_file(),
        mobile_element_references = "mobile-element-genbanks/{sample}.gbk"
    output:
        breseq_dir = directory("breseq-compare-" + sample_info.get_reference_prefix() + "/data/{sample}"),
        html_dir = directory("breseq-compare-" + sample_info.get_reference_prefix() + "/html/{sample}"),
        done_file = "breseq-compare-" + sample_info.get_reference_prefix() + "/html/{sample}/output.done",
        gd_file = "breseq-compare-" + sample_info.get_reference_prefix() + "/gd/{sample}.gd",
        bam = "breseq-compare-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-compare-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        summary_json = "breseq-compare-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json",
    log: 
        "logs/breseq-predict-mutations-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        gd_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/gd"),
        breseq_options = lambda wildcards: get_config_value('BRESEQ_OPTIONS', DEFAULT_BRESEQ_OPTIONS, wildcards.sample),
    threads: 8
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} {params.breseq_options} -s {input.mobile_element_references} -r {input.main_reference} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """
