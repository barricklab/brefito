try: sample_info
except NameError: 
    include: "load-sample-info.smk"

try: DOWNLOAD_DATA_INCLUDED
except NameError: 
    include: "download-data.smk"
    
# Picks the rule order, so it can't vary between samples
NANOPORE_TRIMMING = str(get_uniform_config_value('NANOPORE_TRIMMING', "none"))
    
if NANOPORE_TRIMMING.upper() == "NONE" or NANOPORE_TRIMMING.upper() == "FALSE" or NANOPORE_TRIMMING.upper() == "F" or NANOPORE_TRIMMING == "0":
    print("Nanopore reads are not trimmed. (DEFAULT)")
    ruleorder: trim_nanopore_reads_no_trim > trim_nanopore_reads_porechop_abi > trim_nanopore_reads_porechop
elif NANOPORE_TRIMMING.upper() == "PORECHOP-ABI":
    print("porechop_ABI used for trimming nanopore reads.")
    ruleorder: trim_nanopore_reads_porechop_abi > trim_nanopore_reads_porechop > trim_nanopore_reads_no_trim
elif NANOPORE_TRIMMING.upper() == "PORECHOP":
    print("porechop used for trimming nanopore reads.")
    ruleorder: trim_nanopore_reads_porechop > trim_nanopore_reads_porechop_abi > trim_nanopore_reads_no_trim
else:
    print("Unknown trimming method requested for --config NANOPORE_TRIMMING.")
    print("Valid options are: 'none' (default), 'porechop', 'porechop-ABI'.")
    exit(0)

# Long reads of every type end up in <type>-reads-trimmed/. NANOPORE_TRIMMING picks how nanopore
# reads of every accuracy (nanopore, nanopore-fast, nanopore-hac, nanopore-sup) are trimmed
# (porechop and porechop_ABI look for ONT adapters only). PacBio reads are not trimmed: HiFi reads come
# adapter-free, and CLR subreads have already had their SMRTbell adapters removed by the
# instrument software.
def get_all_trimmed_nanopore_read_names():
    long_read_files = []
    for sample in sample_info.get_sample_list():
        t = sample_info.get_long_read_type(sample)
        if t:
            long_read_files = long_read_files + [ t + "-reads-trimmed/" + d for d in sample_info.get_long_read_list(sample)]
    return long_read_files

rule trim_pacbio_reads_no_trim:
    input:
        "{pacbio_type}-reads/{sample}.fastq.gz"
    output:
        "{pacbio_type}-reads-trimmed/{sample}.fastq.gz"
    wildcard_constraints:
        pacbio_type=PACBIO_TYPE_REGEX
    threads: 1
    shell:
        "ln -s ../{input} {output}"

rule all_trim_nanopore_reads:
    input:
        get_all_trimmed_nanopore_read_names()
    default_target: True

rule trim_nanopore_reads_no_trim:
    input:
        "{ont_type}-reads/{sample}.fastq.gz"
    output:
        "{ont_type}-reads-trimmed/{sample}.fastq.gz"
    wildcard_constraints:
        ont_type=NANOPORE_TYPE_REGEX
    threads: 1
    shell:
        "ln -s ../{input} {output}"

rule trim_nanopore_reads_porechop:
    input:
        "{ont_type}-reads/{sample}.fastq.gz"
    output:
        "{ont_type}-reads-trimmed/{sample}.fastq.gz"
    wildcard_constraints:
        ont_type=NANOPORE_TYPE_REGEX
    log:
        "logs/porechop-{ont_type}-{sample}.log"
    conda:
        "../envs/porechop.yml"
    threads: 12
    shell:
        "porechop --discard_middle -i {input} -o {output} > {log} 2>&1"


rule trim_nanopore_reads_porechop_abi:
    input:
        "{ont_type}-reads/{sample}.fastq.gz"
    output:
        "{ont_type}-reads-trimmed/{sample}.fastq.gz"
    wildcard_constraints:
        ont_type=NANOPORE_TYPE_REGEX
    log:
        "logs/porechop-{ont_type}-{sample}.log"
    conda:
        "../envs/porechop_abi.yml"
    threads: 12
    shell:
        "porechop_abi --ab_initio --discard_middle -i {input} -o {output} > {log} 2>&1"