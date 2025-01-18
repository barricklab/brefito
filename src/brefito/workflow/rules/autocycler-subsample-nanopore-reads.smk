include: "filter-nanopore-reads.smk"

READ_SUBSET_IDS=["01","02","03","04"] # autocycler subsample creates 4 subsampled fastqs for a given fastq file

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
    shell:
        "autocycler subsample --genome_size {config[genome_size]} --reads {input} --out_dir autocycler-nanopore-reads-subsampled/{wildcards.sample} > {log} 2>&1"
