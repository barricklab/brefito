try: sample_info
except NameError:
    include: "load-sample-info.smk"

# Rewrites each GenomeDiff file so that it keeps only its mutation entries, with a dot (.) in
# the id and parent-ids columns. breseq numbers mutations 1..n and lists the evidence ids for
# each one, so inserting or removing a single mutation renumbers everything after it and a
# version control diff covers the whole file. Cleaned files change only where the mutations do.
#
# Needs a breseq build that has the gdtools REMOVE --clear-ids option. Use --breseq-prerelease,
# or --breseq-path to point at a local build.

rule all_clean_genome_diffs_gdtools:
    input:
        ["genome-diffs-clean/" + s + ".gd" for s in sample_info.get_sample_list()]
    default_target: True

rule _clean_genome_diffs_gdtools:
    input:
        genomediff = "genome-diffs/{sample}.gd"
    output:
        "genome-diffs-clean/{sample}.gd"
    conda:
        BRESEQ_ENV
    log:
        "logs/clean-genome-diffs-gdtools-{sample}.log"
    threads: 1
    shell:
        """
        gdtools REMOVE --clear-ids -o {output} {input.genomediff} > {log} 2>&1
        """
