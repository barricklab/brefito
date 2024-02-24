try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "annotate-genomes.smk"

sample_info.set_reference_prefix("annotated-references")

BRESEQ_OPTIONS = "--max-displayed-reads 40 --junction-score-cutoff 0 --junction-minimum-pos-hash-score 10 --brief-html-output"
if 'BRESEQ_OPTIONS' in brefito_config.keys():
    BRESEQ_OPTIONS = brefito_config['BRESEQ_OPTIONS'] 


rule convert_references_to_genbank:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        genbank = "annotated-" + sample_info.get_reference_prefix() + "-genbank/{sample}.gbk"
    log: 
        "logs/breseq-convert-references-to-genbank-{sample}.log"
    conda:
        "../envs/breseq.yml"
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
        "../envs/breseq.yml"
    params:
        reads = "simulated-reads-" + sample_info.get_reference_prefix() + "/{sample}.fastq"
    threads: 8
    shell:
        """
        breseq simulate-reads -m tiled -l 200 -c 80 -r {input.references} -o {params.reads}  > {log} 2>&1
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
        "../envs/breseq.yml"
    params:
        gd_dir = directory("breseq-" + sample_info.get_reference_prefix() + "/gd"),
    threads: 8
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} {BRESEQ_OPTIONS} -s {input.mobile_element_references} -r {input.main_reference} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """
