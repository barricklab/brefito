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

BRACKEN_READ_LENGTH = "100"
if 'BRACKEN_READ_LENGTH' in brefito_config.keys():
    BRACKEN_READ_LENGTH = brefito_config['BRACKEN_READ_LENGTH']


#Default -t, threshold option is 10, this makes it explicit
BRACKEN_OPTIONS = "-r 100 -l S -t 10"
BRACKEN_OPTIONS = "S"
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
        "k2 clean --stop-daemon"

rule classify_kraken2:
    input:
        lambda wildcards: ["illumina-reads-trimmed/" + r for r in sample_info.get_illumina_read_list(wildcards.sample)]
    output:
        report = "classify-kraken2/kracken2_report/{sample}.txt"
    log:
        "logs/classify-kraken2-{sample}.log"
    conda:
        "../envs/kraken2.yml"
    threads: 6
    resources:
        singleton=1
    shell:
        "k2 classify --use-daemon --db {KRAKEN2_DB} --threads {threads} --report {output.report} --output /dev/null  {input} > {log} 2>&1"

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
        "bracken -d {KRAKEN2_DB} -i {input} -o {output} -l {wildcards.classification_level} {BRACKEN_OPTIONS} > {log} 2>&1"

