try: sample_info
except NameError: 
    include: "load-sample-info.smk"

# NOTE: deliberately do NOT include predict-mutations-breseq.smk. This is a
# read-only view over an existing breseq run; if predict were a producer in the
# DAG, a stale-looking breseq output would re-run breseq (often with no reads),
# fail, and Snakemake would delete breseq-<ref>/data/<clone>/. BRESEQ_ENV comes
# from load-sample-info.smk. Missing breseq output is a clean error, not a rebuild.

print()
print("Config options for workflow :: bam2cov-breseq")

BAM2COV_REGION = ""
if 'BAM2COV_REGION' in brefito_config.keys():
    BAM2COV_REGION = brefito_config['BAM2COV_REGION'] 
elif 'REGION' in brefito_config.keys():
    BAM2COV_REGION = brefito_config['REGION'] 
else:
    print()
    print("  You must provide a region for bam2cov. Add to your command: --config BAM2COV_REGION=seq_id:start-end")
    exit(1)


BAM2COV_FORMAT = "png"
if 'BAM2COV_FORMAT' in brefito_config.keys() and brefito_config['BAM2COV_FORMAT'] != None:
    BAM2COV_FORMAT = brefito_config['BAM2COV_FORMAT']
else:
    print()
    print("  Using default --config BAM2COV_FORMAT setting: --format " + BAM2COV_FORMAT)
    print("  You can change this using --config BAM2COV_FORMAT=\"png|svg|pdf|tsv|csv\"")
BAM2COV_FORMAT = BAM2COV_FORMAT.lower()

BAM2COV_OPTIONS = ""
if 'BAM2COV_OPTIONS' in brefito_config.keys() and brefito_config['BAM2COV_OPTIONS'] != None:
    BAM2COV_OPTIONS = brefito_config['BAM2COV_OPTIONS']
    pattern = r"--format +\S+"
    if re.search(pattern, BAM2COV_OPTIONS):
       print()
       print("Provide --format options through --config BAM2COV_FORMAT=\"png|svg|pdf|tsv|csv\" instead of BAM2COV_OPTIONS")
       exit(1)
else:
    print()
    print("  Using default --config BAM2COV_OPTIONS options: " + BAM2COV_OPTIONS)
    print("  You can change this using --config BAM2COV_OPTIONS=\"<options>\"")

print()


BAM2COV_OUTPUT_FILENAME = BAM2COV_REGION + "." + BAM2COV_FORMAT 

rule all_bam2cov_breseq:
    input:
        ["bam2cov-" + sample_info.get_reference_prefix() + "/" + s + "/" + BAM2COV_OUTPUT_FILENAME for s in sample_info.get_sample_list()]
    default_target: True

rule bam2cov_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
    output:
        "bam2cov-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2COV_OUTPUT_FILENAME
    params:
        region = "bam2cov-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2COV_REGION
    log: 
        "logs/bam2cov-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2COV_REGION + ".log"
    conda:
        BRESEQ_ENV
    threads: 1 
    shell:
        """
        echo "breseq BAM2COV --output {output} --fasta {input.fasta} --bam {input.bam} --format {BAM2COV_FORMAT} {BAM2COV_OPTIONS} --region {BAM2COV_REGION}" > {log}
        breseq BAM2COV --output {output} --fasta {input.fasta} --bam {input.bam} --format {BAM2COV_FORMAT} {BAM2COV_OPTIONS} --region {BAM2COV_REGION} >> {log} 2>&1
        """
