try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "predict-mutations-breseq.smk"

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

print("Config options for workflow :: tabulate-ssrs-breseq-*")
print()
SSR_MINIMUM_LENGTH = 6
if 'SSR_MINIMUM_LENGTH' in brefito_config.keys():
    SSR_MINIMUM_LENGTH = brefito_config['SSR_MINIMUM_LENGTH'] 
    print("  User set --minimum-length to " + str(SSR_MINIMUM_LENGTH))
else:
    print("  Using default --minimum-length " + str(SSR_MINIMUM_LENGTH))
    print("  (You can change this using --config SSR_MINIMUM_LENGTH=<int>")

print()

SSR_STRICT_MODE = True
if 'SSR_STRICT_MODE' in brefito_config.keys():
    SSR_STRICT_MODE = str2bool(str(brefito_config['SSR_STRICT_MODE']))
    print("  User set --strict mode to " + str(SSR_STRICT_MODE))
else:
    print("  Using default --strict mode of True")
    print("  You can change this using --config SSR_STRICT_MODE=False")

SSR_STRICT_MODE_ARG = ""
if (SSR_STRICT_MODE):
    SSR_STRICT_MODE_ARG= "--strict"

print()

import os.path

rule all_tabulate_ssrs_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/ssrs/" + s + ".csv" for s in sample_info.get_sample_list()]
    default_target: True

rule tabulate_ssrs_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "breseq-" + sample_info.get_reference_prefix() + "/ssrs/{sample}.csv"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r ')
    log: 
        "logs/breseq-" + sample_info.get_reference_prefix() + "-{sample}-tabulate-ssrs.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq CL-TABULATE {SSR_STRICT_MODE_ARG} --minimum-length {SSR_MINIMUM_LENGTH} -o {output} --bam {input.bam} --fasta {input.fasta} {params.reference_arguments} > {log} 2>&1
        """
