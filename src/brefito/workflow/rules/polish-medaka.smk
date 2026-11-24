include: "filter-nanopore-reads.smk"

# medaka's models are trained on nanopore reads only, so it cannot polish with PacBio reads.
# Stop with a clear message rather than running it on reads it was not built for.
pacbio_samples = [s for s in sample_info.get_sample_list()
                  if sample_info.get_long_read_type(s) not in (None, "nanopore")]
if pacbio_samples:
    print("ERROR: polish-medaka polishes with nanopore reads only, but these samples have PacBio reads:")
    print("  " + " ".join(sorted(pacbio_samples)))
    print("  Polish them with Illumina reads instead (polish-polypolish, polish-polca or polish-breseq).")
    exit(1)

rule polish_with_medaka:
    input:
        reference = "assemblies/{sample}.fasta",
        reads = "nanopore-reads-filtered/{sample}.fastq"
    output:
        path = temp(directory("medaka-polish/{sample}")),
        fasta = "assemblies/{sample}.fasta.polished"
    log: 
        "logs/polish_medaka_{sample}.log"
    conda:
        "../envs/medaka.yml"
    threads: 8
    shell:
      """
      # Copy fasta to output directory in order to avoid preserving an index file
      mkdir -p {output.path}
      cp {input.reference} {output.path}
      medaka_consensus -i {input.reads} -d {output.path}/{wildcards.sample}.fasta -o {output.path} -m r941_min_sup_g507 -t {threads} > {log} 2>&1
      mv {output.path}/consensus.fasta {output.fasta}
      """    