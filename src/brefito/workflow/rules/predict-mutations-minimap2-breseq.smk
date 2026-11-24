try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "align-reads.smk"

import os.path

BRESEQ_THREADS = int(get_uniform_config_value('BRESEQ_THREADS', 8))

# BRESEQ_OPTIONS and NO_DEFAULT_BRESEQ_OPTIONS can be set per sample in data.csv
def get_breseq_options(sample):
    return get_config_value('BRESEQ_OPTIONS', "", sample)

# The sample's long reads (nanopore, pacbio-clr or pacbio-hifi), aligned by minimap2 with the
# preset for their type (align-reads.smk)
def find_available_read_files(wildcards):
    long_read_type = sample_info.get_long_read_type(wildcards.sample)
    if long_read_type is None:
        return []
    prefix = long_read_bam_prefix(long_read_type)
    return [ "aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + wildcards.sample + "/" + prefix + "." + d + ".sam" for d in sample_info.get_long_read_base_list(wildcards.sample)]

# -x (--nanopore) tunes breseq for the homopolymer indels of nanopore reads, so it is only the
# default for nanopore reads. PacBio reads get --require-match-fraction 0.0 alone, since minimap2
# has already decided how much of each read aligns.
def get_breseq_args(sample):
    if config_is_true("NO_DEFAULT_BRESEQ_OPTIONS", sample=sample):
        return ""
    if sample_info.get_long_read_type(sample) == "nanopore":
        return "-x --require-match-fraction 0.0"
    return "--require-match-fraction 0.0"

rule all_predict_mutations_minimap2_breseq:
    input:
        ["minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/" + s + "/output.done" for s in sample_info.get_samples_with_long_reads()]
    default_target: True

rule predict_mutations_minimap2_breseq:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards),
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        breseq_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}"),
        html_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/{sample}"),
        done_file = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/html/{sample}/output.done",
        gd_file = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/gd/{sample}.gd",
        bam = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        summary_json = "minimap2-breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json",
    log: 
        "logs/minimap2-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        gd_dir = directory("minimap2-breseq-" + sample_info.get_reference_prefix() + "/gd"),
        automatic_breseq_args = lambda wildcards: get_breseq_args(wildcards.sample),
        breseq_options = lambda wildcards: get_breseq_options(wildcards.sample),
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    threads: BRESEQ_THREADS
    shell:
        """
        # Create outer directories for moved files
        mkdir -p {params.gd_dir}

        breseq -j {threads} --aligned-sam {params.automatic_breseq_args} {params.breseq_options} {params.reference_arguments} -o {output.breseq_dir}  {input.reads} > {log} 2>&1
        
        # Copy/move output files
        cp {output.breseq_dir}/output/output.gd {output.gd_file}
        rm -rf {output.html_dir}
        mv {output.breseq_dir}/output {output.html_dir}
        """