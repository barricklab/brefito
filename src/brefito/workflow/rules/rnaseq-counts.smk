try: sample_info
except NameError: 
    include: "load-sample-info.smk"

brefito_config['SORT'] = 0
brefito_config['MERGE'] = 1
brefito_config['INDEX'] = 0

HTSEQ_TYPE = "CDS"
if 'HTSEQ_TYPE' in brefito_config.keys():
    HTSEQ_TYPE = brefito_config['HTSEQ_TYPE']
    print("User specified htseq feature type to count " + HTSEQ_TYPE)
else:
    print("Used default htseq feature type to count " + HTSEQ_TYPE)
    print("To change, set --config HTSEQ_TYPE=gene of similar")

# We want bowtiew2 to report one mapping for multimapped pairs
brefito_config['BOWTIE2_OPTIONS'] = "-k 1"

include: "align-reads.smk"

rule merge_rnaseq_counts:
    input:
        ["rnaseq-counts/samples/" + s + ".tsv" for s in sample_info.get_sample_list()]
    output:
        "rnaseq-counts/merged.tsv"
    default_target: True
    run:
        locus_tags = []
        locus_tags_loaded = False
        sample_names = []
        counts = []
        i=0
        input_filename_list = sorted(input) 
        for file_index in range(len(input_filename_list)):
            filename = input_filename_list[file_index]
            print(filename)
            with open(filename) as f:
                sample_name = filename.split(".")[0].split("/")[-1]
                sample_names.append(sample_name)
                locus_tag_index = 0
                for line in f:
                    locus_tag, count = line.rstrip("\n").split("\t")
                    if not locus_tags_loaded:
                        counts.append([])
                    counts[locus_tag_index].append(count)
                    if locus_tags_loaded:
                        if locus_tags[locus_tag_index] != locus_tag:
                            raise ValueError(f"Sample {sample_name}: unexpected locus_tag {locus_tag}")
                    else:
                        locus_tags.append(locus_tag)
                    locus_tag_index += 1
            locus_tags_loaded = True

        with open(output[0], "w") as out:
            out.write(f"locus_tags\t" + "\t".join(sample_names) + "\n")
            for locus_tag_index in range(len(locus_tags)):
                out.write(f"{locus_tags[locus_tag_index]}\t" + "\t".join(counts[locus_tag_index]) + "\n")

rule htseq_rnaseq_counts:
    input:
        gff3 = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.gff3",
        bam =  "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads_merged.bam"
    output: 
        "rnaseq-counts/samples/{sample}.tsv" 
    log:
        "logs/htseq-rnaseq-counts/{sample}.log"
    conda:
        "../envs/htseq.yml"
    shell:
        """
        htseq-count --nonunique random --type {HTSEQ_TYPE} -i ID -r name {input.bam} {input.gff3} > {output}
        """
