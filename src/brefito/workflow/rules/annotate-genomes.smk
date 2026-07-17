# There are some bugs generating Genbank files
# with versions of breseq < 0.38.3, 
# which is why this makes GFF3 only right now!

try: sample_info
except NameError:
    include: "load-sample-info.smk"

# --config SKIP_ISESCAN[=1] (or NO_ISESCAN) skips ISEScan integration; the final .gbk path is unchanged.
# --config SKIP_PHISPY[=1] (or NO_PHISPY) skips PhiSpy; the annotated GenBank is copied through unchanged.
SKIP_ISESCAN = config_is_true("SKIP_ISESCAN") or config_is_true("NO_ISESCAN")
SKIP_PHISPY = config_is_true("SKIP_PHISPY") or config_is_true("NO_PHISPY")

rule all_annotate_genomes:
    input:
        ["annotated-" + sample_info.get_reference_prefix() + "/" + s + ".gbk" for s in sample_info.get_sample_list()]
    default_target: True

#Convert the references to fasta format (so we can re-annotate if needed)
rule convert_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = temp("annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "_{sample}_convert_references.log"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample),
    conda:
        BRESEQ_ENV
    shell:
        """
        breseq CONVERT-REFERENCE -o {output.fasta} -f FASTA {params.reference_arguments} > {log} 2>&1
        """

rule annotate_with_prokka:
    input:
        "annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta"
    output:
        dir = temp(directory("annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}")),
        prokka_annotated_reference = temp("annotated-" + sample_info.get_reference_prefix() + "-prokka/{sample}/reference.gff")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-prokka.log"
    conda:
        "../envs/prokka.yml"
    threads: 4
    shell:
        """
        prokka --cpus {threads} --prefix reference --force --outdir {output.dir} {input} > {log} 2>&1
        """

rule annotate_with_isescan:
    input:
        "annotated-" + sample_info.get_reference_prefix() + "-fasta/{sample}.fasta"
    output:
        # ISEScan exits 0 but writes no *.fasta.csv when it finds zero IS elements, so
        # we depend on the (always-created) output directory instead of that csv. The
        # combine rule checks at runtime whether the csv exists. mkdir -p is a cheap
        # safety net so the directory-output check never fails.
        dir = temp(directory("annotated-" + sample_info.get_reference_prefix() + "-isescan/{sample}"))
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-isescan.log"
    conda:
        "../envs/isescan.yml"
    threads: 8
    shell:
        """
        isescan.py --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1
        mkdir -p {output.dir}
        """

#       isescan.py --removeShortIS --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1


def combine_annotation_inputs(wildcards):
    prefix = "annotated-" + sample_info.get_reference_prefix()
    inputs = {
        "prokka": prefix + "-prokka/" + wildcards.sample + "/reference.gff",
        "prokka_dir": prefix + "-prokka/" + wildcards.sample,
    }
    # Only depend on (and run) ISEScan when it is not being skipped.
    if not SKIP_ISESCAN:
        inputs["isescan_dir"] = prefix + "-isescan/" + wildcards.sample
    return inputs

rule combine_annotation_with_breseq:
    input:
        unpack(combine_annotation_inputs)
    output:
        # Intermediate Prokka(+ISEScan) GenBank; PhiSpy consumes it and produces the
        # final annotated-<prefix>/{sample}.gbk. temp() is safe because PhiSpy copies it.
        temp("annotated-" + sample_info.get_reference_prefix() + "-prokka-isescan/{sample}.gbk")
    params:
        # Expected ISEScan csv path (empty when SKIP_ISESCAN). The shell adds -s only when
        # this file actually exists and is non-empty, so a zero-IS-element result (no csv)
        # falls back to a Prokka-only conversion instead of crashing.
        isescan_csv = lambda wildcards: "" if SKIP_ISESCAN else (
            "annotated-" + sample_info.get_reference_prefix() + "-isescan/" + wildcards.sample + "/" +
            "annotated-" + sample_info.get_reference_prefix() + "-fasta/" + wildcards.sample + ".fasta.csv")
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-combine-annotation-with-breseq.log"
    conda:
        BRESEQ_ENV
    shell:
        """
        if [ -n "{params.isescan_csv}" ] && [ -s "{params.isescan_csv}" ]; then
            breseq CONVERT-REFERENCE -f GENBANK -s "{params.isescan_csv}" -o {output} {input.prokka} > {log} 2>&1
        else
            breseq CONVERT-REFERENCE -f GENBANK -o {output} {input.prokka} > {log} 2>&1
        fi
        """

# Run PhiSpy on the annotated GenBank to add prophage annotations. PhiSpy must run
# after Prokka since it needs an annotated GenBank as input. When SKIP_PHISPY is set,
# the rule copies the input GenBank through unchanged instead of running PhiSpy.
rule annotate_with_phispy:
    input:
        gbk = "annotated-" + sample_info.get_reference_prefix() + "-prokka-isescan/{sample}.gbk"
    output:
        gbk = "annotated-" + sample_info.get_reference_prefix() + "/{sample}.gbk",
        dir = temp(directory("annotated-" + sample_info.get_reference_prefix() + "-phispy/{sample}"))
    params:
        skip = "1" if SKIP_PHISPY else ""
    log:
        "logs/" + "annotated-" + sample_info.get_reference_prefix() + "-{sample}-phispy.log"
    conda:
        "../envs/phispy.yml"
    shell:
        """
        if [ -n "{params.skip}" ]; then
            mkdir -p {output.dir}
            cp {input.gbk} {output.gbk}
        else
            PhiSpy.py {input.gbk} -o {output.dir} > {log} 2>&1
            # PhiSpy's GenBank output (default --output_choice 3) is named after the input basename.
            cp {output.dir}/$(basename {input.gbk}) {output.gbk}
        fi
        """
