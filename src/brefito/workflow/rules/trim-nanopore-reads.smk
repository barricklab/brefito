try: sample_info
except NameError: 
    include: "load-sample-info.smk"

try: DOWNLOAD_DATA_INCLUDED
except NameError: 
    include: "download-data.smk"
    
NANOPORE_TRIMMING = "none"
if 'NANOPORE_TRIMMING' in brefito_config.keys():
    NANOPORE_TRIMMING = bool(brefito_config['NANOPORE_TRIMMING'])

if NANOPORE_TRIMMING.upper() == "NONE" or NANOPORE_TRIMMING.upper() == "FALSE" or NANOPORE_TRIMMING.upper() == "F" or NANOPORE_TRIMMING == "0":
    print("Nanopore reads are not trimmed.")
    ruleorder: trim_nanopore_reads_no_trim > trim_nanopore_reads_porechop_abi > trim_nanopore_reads_porechop
elif NANOPORE_TRIMMING.upper() == "PORECHOP-ABI":
    print("porechop_ABI used for trimming nanopore reads.")
    ruleorder: trim_nanopore_reads_porechop_abi > trim_nanopore_reads_porechop_abi > trim_nanopore_reads_no_trim
elif NANOPORE_TRIMMING.upper() == "PORECHOP":
    print("porechop used for trimming nanopore reads. (DEFAULT)")
    ruleorder: trim_nanopore_reads_porechop > trim_nanopore_reads_porechop_abi > trim_nanopore_reads_no_trim
else:
    print("Unknown trimming method requested for --config NANOPORE_TRIMMING.")
    print("Valid options are: 'none', 'porechop', 'porechop-ABI'.")
    exit(0)

def get_all_trimmed_nanopore_read_names():
    nanopore_files = []
    for sample in sample_info.get_sample_list():
        nanopore_files = nanopore_files + [ "nanopore-reads-trimmed/" + d for d in sample_info.get_nanopore_read_list(sample)] 
    return nanopore_files

rule all_trim_nanopore_reads:
    input:
        get_all_trimmed_nanopore_read_names()
    default_target: True

rule trim_nanopore_reads_no_trim:
    input:
        "nanopore-reads/{sample}.fastq.gz"
    output:
        "nanopore-reads-trimmed/{sample}.fastq.gz"
    threads: 1
    shell:
        "ln -s ../{input} {output}"

rule trim_nanopore_reads_porechop:
    input:
        "nanopore-reads/{sample}.fastq.gz"
    output:
        "nanopore-reads-trimmed/{sample}.fastq.gz"
    log:
        "logs/porechop-{sample}.log"
    conda:
        "../envs/porechop.yml"
    threads: 12
    shell:
        "porechop --discard_middle -i {input} -o {output} > {log} 2>&1"


rule trim_nanopore_reads_porechop_abi:
    input:
        "nanopore-reads/{sample}.fastq.gz"
    output:
        "nanopore-reads-trimmed/{sample}.fastq.gz"
    log:
        "logs/porechop-{sample}.log"
    conda:
        "../envs/porechop_abi.yml"
    threads: 12
    shell:
        "porechop_abi --ab_initio --discard_middle -i {input} -o {output} > {log} 2>&1"