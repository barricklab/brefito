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

# REGION and FORMAT name the output files, so they can't vary between samples
BAM2ALN_REGION = get_uniform_config_value('BAM2ALN_REGION', get_uniform_config_value('REGION', ""))
if BAM2ALN_REGION == "":
    print()
    print("  You must provide a region for bam2aln. Add to your command: --config BAM2ALN_REGION=seq_id:start-end")
    exit(1)



BAM2ALN_FORMAT = "HTML"
if get_uniform_config_value('BAM2ALN_FORMAT') != None:
    BAM2ALN_FORMAT = get_uniform_config_value('BAM2ALN_FORMAT')
else:
    print()
    print("  Using default --config BAM2ALN_FORMAT setting: --format " + BAM2ALN_FORMAT)
    print("  You can change this using --config BAM2ALN_FORMAT=\"HTML|TXT|JSON\"")

BAM2ALN_OPTIONS = ""

# BAM2ALN_OPTIONS can be set per sample in data.csv
def get_aln_options(sample):
    options = get_config_value('BAM2ALN_OPTIONS', BAM2ALN_OPTIONS, sample)
    if options == None:
        options = ""
    if re.search(r"--format +\S+", options):
        print()
        print("Provide --format options through --config BAM2ALN_FORMAT=\"HTML|TXT|JSON\" instead of BAM2ALN_OPTIONS")
        exit(1)
    return options

for this_sample in sample_info.get_sample_list():
    get_aln_options(this_sample)

if 'BAM2ALN_OPTIONS' not in brefito_config.keys():
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
        aln_options = lambda wildcards: get_aln_options(wildcards.sample),
        region = "bam2aln-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2ALN_OUTPUT_FILENAME
    log: 
        "logs/bam2aln-" + sample_info.get_reference_prefix() + "/{sample}/" + BAM2ALN_REGION + ".log"
    conda:
        BRESEQ_ENV
    threads: 1
    shell:
        """
        echo "breseq BAM2ALN --output {params.region} --fasta {input.fasta} --bam {input.bam} --format {BAM2ALN_FORMAT} {params.aln_options} --region {BAM2ALN_REGION}" > {log} 2>&1
        breseq BAM2ALN --output {params.region} --fasta {input.fasta} --bam {input.bam} --format {BAM2ALN_FORMAT} {params.aln_options} --region {BAM2ALN_REGION} >> {log} 2>&1
        """
