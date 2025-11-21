include: "predict-mutations-breseq.smk"

TILE_SIZE = 10000
if 'TILE_SIZE' in brefito_config.keys():
    TILE_SIZE = int(brefito_config['TILE_SIZE'])

TILE_OVERLAP = ""
if 'TILE_OVERLAP' in brefito_config.keys():
    TILE_OVERLAP = brefito_config['TILE_OVERLAP'] 

PLOT_RESOLUTION = 1000
if 'PLOT_RESOLUTION' in brefito_config.keys():
    PLOT_RESOLUTION = bool(brefito_config['PLOT_RESOLUTION'])

PLOT_FORMAT = "PNG"
if 'PLOT_FORMAT' in brefito_config.keys():
    PLOT_FORMAT = bool(brefito_config['PLOT_FORMAT'])

print("Using breseq bam2cov options, --tile-size {} --tile-overlap {} --resolution {} --format {}".format(TILE_SIZE, TILE_OVERLAP, PLOT_RESOLUTION, PLOT_FORMAT))
print("You can set these using --config TILE_SIZE=<int> --config TILE_OVERLAP=<int> --config PLOT_RESOLUTION=<int> --config PLOT_FORMAT=<PDF/PNG>\n")

rule all_coverage_plots_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/cov/" + s + "/coverage.done" for s in sample_info.get_sample_list()]
    default_target: True

# We use the breseq files to get the average coverage
# Need to update bam2cov to create  output files for each reference by default if no argument provided...
rule coverage_plots_breseq:
    input:
        bam = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.bam",
        fasta = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/reference.fasta",
        breseq_summary_json = "breseq-" + sample_info.get_reference_prefix() + "/data/{sample}/data/summary.json"
    output:
        dir = directory("breseq-" + sample_info.get_reference_prefix() + "/cov/{sample}"),
        #whole_genome_file = "evaluate/coverage_plots/{technology}/{sample}/all.pdf",
        done_file = "breseq-" + sample_info.get_reference_prefix() + "/cov/{sample}/coverage.done"
    log:
        "logs/coverage-plots-breseq-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/breseq.yml"
    threads: 1
    shell:
        """
        breseq bam2cov -o {output.dir} -x tiling_ -b {input.bam} -f {input.fasta} -a --format {PLOT_FORMAT} --tile-size {TILE_SIZE} --tile-overlap {TILE_OVERLAP} --resolution {PLOT_RESOLUTION} -j {input.breseq_summary_json} > {log} 2>&1
        breseq bam2cov -o {output.dir} -x summary_ -b {input.bam} -f {input.fasta} -a -j {input.breseq_summary_json} -p 10000 >> {log} 2>&1
        touch {output.done_file}
        """
