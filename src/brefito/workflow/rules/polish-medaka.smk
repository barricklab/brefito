include: "filter-nanopore-reads.smk"

# medaka's models are trained on nanopore reads only, so it cannot polish with PacBio reads.
# Stop with a clear message rather than running it on reads it was not built for.
pacbio_samples = [s for s in sample_info.get_sample_list()
                  if sample_info.get_long_read_type(s) is not None
                  and not is_nanopore_type(sample_info.get_long_read_type(s))]
if pacbio_samples:
    print("ERROR: polish-medaka polishes with nanopore reads only, but these samples have PacBio reads:")
    print("  " + " ".join(sorted(pacbio_samples)))
    print("  Polish them with Illumina reads instead (polish-polypolish or polish-breseq).")
    exit(1)

# The medaka model must match the flow cell and basecaller of the reads. The default comes from
# the read type (LONG_READ_SETTINGS): an R9.4.1 SUP model for nanopore and nanopore-fast reads,
# an R10.4.1 HAC model for nanopore-hac reads and an R10.4.1 SUP model for nanopore-sup reads.
# MEDAKA_MODEL (per sample in data.csv, or
# --config) sets it explicitly; 'medaka tools list_models' lists the choices.
def get_medaka_model(sample):
    default_model = long_read_setting(sample_info.get_long_read_type(sample), "medaka_model")
    return get_config_value('MEDAKA_MODEL', default_model, sample)

rule polish_with_medaka:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = lambda wildcards: sample_info.get_long_read_path(wildcards.sample, "filtered")
    output:
        path = temp(directory("medaka-polish/{sample}")),
        fasta = "assemblies/{sample}.fasta.polished"
    log:
        "logs/polish_medaka_{sample}.log"
    conda:
        "../envs/medaka.yml"
    params:
        model = lambda wildcards: get_medaka_model(wildcards.sample)
    threads: 8
    shell:
      """
      # Copy fasta to output directory in order to avoid preserving an index file
      mkdir -p {output.path}
      cp {input.reference} {output.path}
      medaka_consensus -i {input.reads} -d {output.path}/{wildcards.sample}.fasta -o {output.path} -m {params.model} -t {threads} > {log} 2>&1
      mv {output.path}/consensus.fasta {output.fasta}
      """
