include: "filter-nanopore-reads.smk"

FLYE_RESAMPLED_IDS = ["01", "07", "13", "19"]
MINIASM_MINIPOLISH_RESAMPLED_IDS = ["02", "08", "14", "20"]
RAVEN_RESAMPLED_IDS = ["03", "09", "15", "21"]
CANU_RESAMPLED_IDS = ["04", "10", "16", "22"]
NECAT_RESAMPLED_IDS = ["05", "11", "17", "23"]
UNICYCLER_RESAMPLED_IDS = ["06", "12", "18", "24"]

ALL_RESAMPLED_IDS = FLYE_RESAMPLED_IDS + RAVEN_RESAMPLED_IDS + MINIASM_MINIPOLISH_RESAMPLED_IDS + CANU_RESAMPLED_IDS + NECAT_RESAMPLED_IDS + UNICYCLER_RESAMPLED_IDS

rule subsample_nanopore_reads:
    input:
        "02_filtered_nanopore_reads/{sample}.fastq"
    output:
        expand("03_subsampled_nanopore_reads/{{sample}}/sample_{assembly_id}.fastq", assembly_id=ALL_RESAMPLED_IDS)
    conda:
        "../envs/trycycler.yml"
    shell:
        "trycycler subsample --genome_size {config[genome_size]} --reads {input} --count 24 --out_dir 03_subsampled_nanopore_reads/{wildcards.sample}"


