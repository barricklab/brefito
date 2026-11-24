try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "filter-nanopore-reads.smk"


# The default is the read-type option for the sample's long reads: --nano-raw (nanopore),
# --pacbio-raw (PacBio CLR) or --pacbio-hifi (PacBio HiFi). FLYE_OPTIONS replaces it, so it must
# include the read-type option itself.
DEFAULT_FLYE_OPTIONS = "--nano-raw"
def get_default_flye_options(sample):
    return long_read_setting(sample_info.get_long_read_type(sample), "flye")

if 'FLYE_OPTIONS' in brefito_config.keys():
    print("flye assembly with user-specified command-line options:")
    print("  " + brefito_config['FLYE_OPTIONS'])
else:
    print("flye assembly with default command-line options (unless set per sample in data.csv):")
    print("  " + DEFAULT_FLYE_OPTIONS + "   (PacBio CLR: --pacbio-raw, PacBio HiFi: --pacbio-hifi)")
    print("To change these options add something like this to your brefito command")
    print("  --config FLYE_OPTIONS=\"--meta --nano-hq\"")

rule all:
    input: [ "assemblies/{}.fasta".format(d) for d in sample_info.get_samples_with_long_reads()]

rule assemble_with_flye:
    input:
        lambda wildcards: sample_info.get_long_read_path(wildcards.sample, "filtered")
    output:
        fasta = "assemblies/{sample}.fasta",
        gfa = "assemblies/{sample}.gfa",
        tmpdir = temp(directory("flye_assembly/{sample}"))
    log:
        "logs/flye_{sample}.log"
    conda:
        "../envs/flye.yml"
    params:
        flye_options = lambda wildcards: get_config_value('FLYE_OPTIONS', get_default_flye_options(wildcards.sample), wildcards.sample)
    threads: 8
    shell:
        """
        mkdir -p {output.tmpdir}
        flye {params.flye_options} {input} --threads {threads} --out-dir {output.tmpdir} > {log} 2>&1
        cp {output.tmpdir}/assembly.fasta {output.fasta}
        cp {output.tmpdir}/assembly_graph.gfa {output.gfa}
        """
