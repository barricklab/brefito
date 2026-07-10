try: sample_info
except NameError:
    include: "load-sample-info.smk"

OUTPUT_FORMAT="GENBANK"

if "OUTPUT_FORMAT" in brefito_config.keys():
    OUTPUT_FORMAT = brefito_config['OUTPUT_FORMAT'].upper()
    print(f"User specified output format is {OUTPUT_FORMAT}")
else:
    print("Using default output format: " + OUTPUT_FORMAT)
    print("You can change this by providing --config OUTPUT_FORMAT=<GFF3|GENBANK|FASTA>")

if not(OUTPUT_FORMAT in ['GFF3', 'GENBANK', 'FASTA']):
    raise ValueError(f"File format {OUTPUT_FORMAT} not recognized. Please provide the output format. Valid args are GFF3, GENBANK, FASTA")

file_format_ext = {'GFF3':'.gff3', 'GENBANK':'.gbk', 'FASTA':'.fasta'}

rule all_mutate_genomes_gdtools:
    input:
        ["mutants/" + s + file_format_ext[OUTPUT_FORMAT] for s in sample_info.get_sample_list()]
    default_target: True

rule mutate_genomes_gdtools:
    input:
        genomediff = "genome-diffs/{sample}.gd",
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        "mutants/{sample}" + file_format_ext[OUTPUT_FORMAT]
    log:
        "logs/mutate-genomes-gdtools-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        BRESEQ_ENV
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample, '-r '),
        output_format = OUTPUT_FORMAT
    threads: 1
    shell:
        """
        gdtools APPLY -o {output} {params.reference_arguments} -f {params.output_format}  {input.genomediff} -r {input.references} > {log} 2>&1
        """
