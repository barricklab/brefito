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

for this_sample in sample_info.get_sample_list():
    get_genome_size(this_sample)

rule all_subsample_nanopore_reads:
    input:
        ["autocycler-nanopore-reads-subsampled/" + s + "/sample_01.fastq" for s in sample_info.get_sample_list() ]
    default_target: True

rule subsample_nanopore_reads:
    input:
        "nanopore-reads-filtered/{sample}.fastq"
    output:
        expand("autocycler-nanopore-reads-subsampled/{{sample}}/sample_{subset_id}.fastq", subset_id=READ_SUBSET_IDS)
    conda:
        "../envs/autocycler.yml"    
    log:
        "logs/autocycler-subsample-{sample}.log"
    params:
        genome_size = lambda wildcards: get_genome_size(wildcards.sample)
    shell:
        "autocycler subsample --genome_size {params.genome_size} --reads {input} --out_dir autocycler-nanopore-reads-subsampled/{wildcards.sample} > {log} 2>&1"
