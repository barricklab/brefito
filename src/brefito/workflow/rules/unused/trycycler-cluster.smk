include: "trycycler-assemble.smk"

rule trycycler_cluster:
    input:
        reads = "01_trimmed_reads/{dataset}.long_reads.trimmed.fastq",
        assemblies = expand("04_assemblies/{{dataset}}/assembly_{assembly_id}.fasta", assembly_id=ALL_RESAMPLED_IDS)
    output:
        output_directory = directory("05_trycycler/{dataset}"),
        done_file = "05_trycycler/{dataset}/clustering.done"
    conda:
        "../envs/trycycler.yml"
    threads: 128
    shell:
        "trycycler cluster --threads {threads} --assemblies {input.assemblies} --reads {input.reads} --out_dir 05_trycycler/{wildcards.dataset} && touch {output.done_file}"

