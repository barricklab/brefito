# This version of classify-kraken2 copies the hash file
# into /dev/shm RAM so it can spawn multiple processes
# that all use the shared memory for quick execution.

import os

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

# Different options

if 'KRAKEN2_DB' in brefito_config.keys():
    KRAKEN2_DB = brefito_config['KRAKEN2_DB']
else:
    raise Exception("Must provide --config KRAKEN2_DB=<path to database>")

# Location in shared memory to copy
SHM_DB = f"/dev/shm/{os.environ['USER']}/k2db"

BRACKEN_READ_LENGTH = "100"
if 'BRACKEN_READ_LENGTH' in brefito_config.keys():
    BRACKEN_READ_LENGTH = brefito_config['BRACKEN_READ_LENGTH']


#Default -t, threshold option is 10, this makes it explicit
BRACKEN_OPTIONS = "-r " + BRACKEN_READ_LENGTH + " -t 10"
if 'BRACKEN_OPTIONS' in brefito_config.keys():
    BRACKEN_OPTIONS = brefito_config['BRACKEN_OPTIONS']

classification_levels=["P","C","O","F","G","S"]


rule all_classify_kraken2:
    input:
        expand("classify-kraken2/bracken_output/{sample}_classify_{classification_level}.txt", sample = sample_info.sample_list, classification_level=classification_levels)
    default_target: True
    conda:
        "../envs/kraken2.yml"
    shell:
        "rm -rf {SHM_DB}"

onerror:
    shell(f"rm -rf {SHM_DB} k2db_in_shm || true")

rule kraken2_db_to_shm:
    input:
        hash=f"{KRAKEN2_DB}/hash.k2d",
        opts=f"{KRAKEN2_DB}/opts.k2d",
        taxo=f"{KRAKEN2_DB}/taxo.k2d",
    output:
        temp(touch("k2db_in_shm"))
    params:
        dest=SHM_DB,
    conda:
        "../envs/kraken2.yml"
    shell:
        r"""
        rm -rf {SHM_DB}
        mkdir -p {SHM_DB}
        cp -f {input} {SHM_DB}/
        #k2 inspect --db {SHM_DB} --skip-counts --memory-mapping > /dev/null
        """

rule classify_kraken2:
    input:
        lambda wildcards: ["illumina-reads-trimmed/" + r for r in sample_info.get_illumina_read_list(wildcards.sample)],
        "k2db_in_shm"
    output:
        report = "classify-kraken2/kracken2_report/{sample}.txt"
    log:
        "logs/classify-kraken2-{sample}.log"
    conda:
        "../envs/kraken2.yml"
    threads: 8
    shell:
        "k2 classify --memory-mapping --db {SHM_DB} --threads {threads} --report {output.report} --output /dev/null  {input} > {log} 2>&1"

rule classify_bracken:
    input:
        "classify-kraken2/kracken2_report/{sample}.txt"
    output:
        "classify-kraken2/bracken_output/{sample}_classify_{classification_level}.txt"
    log:
        "logs/classify-bracken-{sample}-{classification_level}.log"
    conda:
        "../envs/kraken2.yml"
    threads: 1
    shell:
        "bracken -d {KRAKEN2_DB} -i {input} -o {output} {BRACKEN_OPTIONS} -l {wildcards.classification_level} > {log} 2>&1"

