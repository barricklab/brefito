include: "filter-nanopore-reads.smk"

READ_SUBSET_IDS=["01","02","03","04"] # autocycler subsample creates 4 subsampled fastqs for a given fastq file

# GENOME_SIZE is required: --config genome_size=<bp>, or per sample in data.csv
def get_genome_size(sample):
    genome_size = get_config_value('GENOME_SIZE', None, sample)
    if genome_size == None or str(genome_size).strip() == "":
        print("ERROR: No genome size for sample " + sample + ". Add --config genome_size=<bp>")
        print("  or a '" + sample + ",genome_size,<bp>' row to data.csv.")
        exit(1)
    return genome_size

# Only samples with long reads are assembled, so only they need a genome size
for this_sample in sample_info.get_samples_with_long_reads():
    get_genome_size(this_sample)

# Long reads of any type: <type>-reads-filtered/ -> autocycler-<type>-reads-subsampled/
def get_autocycler_subsampled_reads(sample, subset_id):
    return "autocycler-" + sample_info.get_long_read_type(sample) + "-reads-subsampled/" + sample + "/sample_" + subset_id + ".fastq"

rule all_subsample_nanopore_reads:
    input:
        [get_autocycler_subsampled_reads(s, "01") for s in sample_info.get_samples_with_long_reads()]
    default_target: True

rule subsample_nanopore_reads:
    input:
        "{lrt}-reads-filtered/{sample}.fastq"
    output:
        expand("autocycler-{{lrt}}-reads-subsampled/{{sample}}/sample_{subset_id}.fastq", subset_id=READ_SUBSET_IDS)
    wildcard_constraints:
        lrt=LONG_READ_TYPE_REGEX
    conda:
        "../envs/autocycler.yml"
    log:
        "logs/autocycler-subsample-{lrt}-{sample}.log"
    params:
        genome_size = lambda wildcards: get_genome_size(wildcards.sample)
    shell:
        "autocycler subsample --genome_size {params.genome_size} --reads {input} --out_dir autocycler-{wildcards.lrt}-reads-subsampled/{wildcards.sample} > {log} 2>&1"
