include: "filter-nanopore-reads.smk"

FLYE_RESAMPLED_IDS = ["01", "07", "13", "19"]
MINIASM_MINIPOLISH_RESAMPLED_IDS = ["02", "08", "14", "20"]
RAVEN_RESAMPLED_IDS = ["03", "09", "15", "21"]
CANU_RESAMPLED_IDS = ["04", "10", "16", "22"]
NECAT_RESAMPLED_IDS = ["05", "11", "17", "23"]
UNICYCLER_RESAMPLED_IDS = ["06", "12", "18", "24"]

ALL_RESAMPLED_IDS = FLYE_RESAMPLED_IDS + RAVEN_RESAMPLED_IDS + MINIASM_MINIPOLISH_RESAMPLED_IDS + CANU_RESAMPLED_IDS + NECAT_RESAMPLED_IDS + UNICYCLER_RESAMPLED_IDS

SUBSAMPLE_OPTIONS = "--genome_size " + str(config["genome_size"])
if "SUBSAMPLE_OPTIONS" not in config.keys():
    print("Subsampling nanopore reads using default trycycler subsample command-line options:")
    print("  " + SUBSAMPLE_OPTIONS)
    print("To change these options add something like this to your brefito command")
    print("  --config SUBSAMPLE_OPTIONS=\"--genome_size 4600000 --min_read_depth 60\"")
else:
    SUBSAMPLE_OPTIONS = config["SUBSAMPLE_OPTIONS"]
    print("Subsampling nanopore reads with user-specified trycycler subsample command-line options:")
    print("  " + SUBSAMPLE_OPTIONS)

rule all_subsample_nanopore_reads:
    input:
        ["nanopore-reads-subsampled/" + s + "/sample_01.fastq" for s in sample_info.get_sample_list() ]
    default_target: True

rule subsample_nanopore_reads:
    input:
        "nanopore-reads-filtered/{sample}.fastq"
    output:
        expand("nanopore-reads-subsampled/{{sample}}/sample_{assembly_id}.fastq", assembly_id=ALL_RESAMPLED_IDS)
    conda:
        "../envs/trycycler.yml"
    log:
        "logs/trycycler-subsample-{sample}.log"
    shell:
        "trycycler subsample {SUBSAMPLE_OPTIONS} --reads {input} --count 24 --out_dir nanopore-reads-subsampled/{wildcards.sample} > {log} 2>&1"


