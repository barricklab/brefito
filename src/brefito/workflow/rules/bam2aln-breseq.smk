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

BAM2ALN_REGION = ""
if 'BAM2ALN_REGION' in brefito_config.keys():
    BAM2ALN_REGION = brefito_config['BAM2ALN_REGION'] 
elif 'REGION' in brefito_config.keys():
    BAM2ALN_REGION = brefito_config['REGION'] 
else:
    print()
    print("  You must provide a region for bam2aln. Add to your command: --config BAM2ALN_REGION=seq_id:start-end")
    exit(1)



BAM2ALN_FORMAT = "HTML"
if 'BAM2COV_FORMAT' in brefito_config.keys():
    BAM2COV_FORMAT = brefito_config['BAM2COV_FORMAT']
else:
    print()
    print("  Using default --config BAM2ALN_FORMAT setting: --format " + BAM2ALN_FORMAT)
    print("  You can change this using --config BAM2ALN_FORMAT=\"HTML|TXT|JSON\"")

BAM2ALN_OPTIONS = ""
if 'BAM2ALN_OPTIONS' in brefito_config.keys():
    BAM2ALN_OPTIONS = brefito_config['BAM2ALN_OPTIONS']
    pattern = r"--format +\S+"
    if re.search(pattern, my_string):
       print()
       print("Provide --format options through --config BAM2ALN_FORMAT=\"HTML|TXT|JSON\" instead of BAM2ALN_OPTIONS")
       exit(1)
else:
    print()
    print("  Using default --config BAM2ALN_OPTIONS options: " + BAM2ALN_OPTIONS)
    print("  You can change this using --config BAM2ALN_OPTIONS=\"<options>\"")

print()


BAM2ALN_OUTPUT_FILENAME = BAM2ALN_REGION + "." + BAM2ALN_FORMAT.lower()

rule all_bam2aln_breseq:
    input:
        ["bam2aln-" + sample_info.get_reference_prefix() + "/" + s + "/" + BAM2ALN_OUTPUT_FILENAME for s in sample_info.get_sample_list()]
    default_target: True

rule bam2aln_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
    output:
        "bam2aln-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2ALN_OUTPUT_FILENAME
    params:
        region = "bam2aln-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2ALN_OUTPUT_FILENAME
    log: 
        "logs/bam2aln-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2ALN_REGION + ".log"
    conda:
        BRESEQ_ENV
    threads: 1
    shell:
        """
        echo "breseq BAM2ALN --output {params.region} --fasta {input.fasta} --bam {input.bam} --format {BAM2ALN_FORMAT} {BAM2ALN_OPTIONS} --region {BAM2ALN_REGION}" > {log} 2>&1
        breseq BAM2ALN --output {params.region} --fasta {input.fasta} --bam {input.bam} --format {BAM2ALN_FORMAT} {BAM2ALN_OPTIONS} --region {BAM2ALN_REGION} >> {log} 2>&1
        """
