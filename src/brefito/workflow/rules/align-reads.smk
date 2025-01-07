try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-nanopore-reads.smk"
include: "trim-illumina-reads.smk"

# Different options
MERGE = 0      # Merge together all of the BAM files
SORT = 1      # Sort reads by coordinate in the BAM file
INDEX = 1     # Create BAM indexes *.bai

if 'SORT' in brefito_config.keys():
    SORT = brefito_config['SORT']

if 'MERGE' in brefito_config.keys():
    MERGE = brefito_config['MERGE']

if 'INDEX' in brefito_config.keys():
    INDEX = brefito_config['INDEX']

# Can't index if we don't sort
if SORT == 0:
    print("Can't INDEX unsorted BAM files. Setting SORT to 0 (false).")
    INDEX=0

#Note: It's too hard to make temp() work when sorting and merging are conditional, so we clean up all prior files within rules

#print("Merge: " + str(MERGE))
#print("Sort: " + str(SORT))
#print("Index: " + str(INDEX))

## Need to fix rule priorities for different sorting and merging cases
ruleorder: samtools_view > samtools_sort

## Change output files 
if MERGE:
    if SORT:
        ## Steps used: samtools_view_for_sort > samtools_sort_for_merge > samtools_merge

        SAMTOOLS_VIEW_INPUT_SUFFIX  = ".sam"
        SAMTOOLS_VIEW_OUTPUT_SUFFIX = ".unmerged.unsorted.bam"

        SAMTOOLS_SORT_INPUT_SUFFIX  = SAMTOOLS_VIEW_OUTPUT_SUFFIX
        SAMTOOLS_SORT_OUTPUT_SUFFIX = ".unmerged.bam"

        SAMTOOLS_MERGE_INPUT_SUFFIX  = SAMTOOLS_SORT_OUTPUT_SUFFIX
        SAMTOOLS_MERGE_OUTPUT_SUFFIX = ".bam"

    else:
        ## Steps used: samtools_view > samtools_merge

        SAMTOOLS_VIEW_INPUT_SUFFIX  = ".sam"
        SAMTOOLS_VIEW_OUTPUT_SUFFIX = ".unmerged.bam"

        SAMTOOLS_SORT_INPUT_SUFFIX  = ".unused.1"
        SAMTOOLS_SORT_OUTPUT_SUFFIX = ".unused.2"

        SAMTOOLS_MERGE_INPUT_SUFFIX  = SAMTOOLS_VIEW_OUTPUT_SUFFIX
        SAMTOOLS_MERGE_OUTPUT_SUFFIX = ".bam"

else: # Not merged

    if SORT:
        ## Steps used: samtools_view_for_sort > samtools_sort
        SAMTOOLS_VIEW_INPUT_SUFFIX  = ".sam"
        SAMTOOLS_VIEW_OUTPUT_SUFFIX = ".unsorted.bam"

        SAMTOOLS_SORT_INPUT_SUFFIX = SAMTOOLS_VIEW_OUTPUT_SUFFIX
        SAMTOOLS_SORT_OUTPUT_SUFFIX =  ".bam"

        SAMTOOLS_MERGE_INPUT_SUFFIX = ".unused.1"
        SAMTOOLS_MERGE_OUTPUT_SUFFIX =  ".unused.2"

    else:
        ## Steps used: samtools_view
        ruleorder: samtools_view > samtools_view_with_sort

        SAMTOOLS_VIEW_INPUT_SUFFIX  = ".sam"
        SAMTOOLS_VIEW_OUTPUT_SUFFIX = ".bam"

        SAMTOOLS_SORT_INPUT_SUFFIX =  ".unused.1"
        SAMTOOLS_SORT_OUTPUT_SUFFIX =   ".unused.2"

        SAMTOOLS_MERGE_INPUT_SUFFIX =  ".unused.3"
        SAMTOOLS_MERGE_OUTPUT_SUFFIX =   ".unused.4"


BOWTIE2_OPTIONS = ""
if 'BOWTIE2_OPTIONS' in brefito_config.keys():
    BOWTIE2_OPTIONS = brefito_config['BOWTIE2_OPTIONS'] 

BOWTIE2_PE_OPTIONS = ""
if 'BOWTIE2_PE_OPTIONS' in brefito_config.keys():
    BOWTIE2_PE_OPTIONS = brefito_config['BOWTIE2_PE_OPTIONS'] 

BOWTIE2_SE_OPTIONS = ""
if 'BOWTIE2_SE_OPTIONS' in brefito_config.keys():
    BOWTIE2_SE_OPTIONS = brefito_config['BOWTIE2_SE_OPTIONS'] 

## The output files change depending on the options

# If not merged, then they are in a folder and named according to

def get_all_output_files_list():
    all_output_files_list = []

    for s in sample_info.get_sample_list():

        all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.fasta")
        if INDEX:
            all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.fasta.fai")
        all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/reference.gff3")
        
        if not MERGE:
            for n in sample_info.get_nanopore_read_base_list(s):
                all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads.{}.bam".format(n))
                if INDEX:
                    all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads.{}.bam.bai".format(n))

            for i_se in sample_info.get_illumina_SE_read_base_list(s):
                all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.SE.bam".format(i_se))
                if INDEX:
                    all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.SE.bam.bai".format(i_se))

            for i_pe in sample_info.get_illumina_PE_read_base_list(s):
                if INDEX:
                    all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.PE.bam".format(i_pe))
                all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads.{}.PE.bam.bai".format(i_pe))

        else: # MERGE
            if len(sample_info.get_nanopore_read_base_list(s)) > 0:
                all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads_merged.bam")
                if INDEX:
                    all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/nanopore_reads_merged.bam.bai")
            if len(sample_info.get_illumina_SE_read_base_list(s) + sample_info.get_illumina_PE_read_base_list(s)) > 0:
                all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads_merged.bam")
                if INDEX:
                    all_output_files_list.append("aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + s + "/illumina_reads_merged.bam.bai")

    #print(all_output_files_list)
    return all_output_files_list

rule all_align_reads:
    input:
        get_all_output_files_list()
    default_target: True

rule align_nanopore_reads:
    input:
        reads = "nanopore-reads-trimmed/{reads}.fastq.gz",
        reference = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta"
    output:
        temp("aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/nanopore_reads.{reads}.sam")
    log:
        "logs/align-nanopore-reads-" + sample_info.get_reference_prefix() + "-{sample}-{reads}.log"
    conda:
        "../envs/minimap2.yml"
    threads: 12
    shell:
        "minimap2 --secondary=no -t {threads} -x map-ont -a -o {output} {input.reference} {input.reads} > {log} 2>&1"

rule align_PE_illumina_reads:
    input:
        reads = expand("illumina-reads-trimmed/{{reads}}.{read_num}.fastq.gz", read_num=["R1", "R2"]),
        reference = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta"
    output:
        temp("aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.PE.sam")
    log:
        "logs/align-illumina-reads-" + sample_info.get_reference_prefix() + "-{sample}-{reads}.log"
    conda:
        "../envs/bowtie2.yml"
    params:
        bowtie2_index = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.PE.bowtie2_index"
    threads: 12
    shell:
        """
        bowtie2-build --threads {threads} {input.reference} {params.bowtie2_index} > {log} 2>&1
        bowtie2 {BOWTIE2_PE_OPTIONS} {BOWTIE2_OPTIONS} --threads {threads} -x {params.bowtie2_index} -1 {input.reads[0]} -2 {input.reads[1]} -S {output} >> {log} 2>&1
        rm {params.bowtie2_index}*
        """

rule align_SE_illumina_reads:
    input:
        reads = expand("illumina-reads-trimmed/{{reads}}.{read_num}.fastq.gz", read_num=["SE"]),
        reference = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta"
    output:
        temp("aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.SE.sam")
    log:
        "logs/align-SE-illumina-reads-" + sample_info.get_reference_prefix() + "-{sample}-{reads}.log"
    conda:
        "../envs/bowtie2.yml"
    params:
        bowtie2_index = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads.{reads}.SE.bowtie2_index"
    threads: 12
    shell:
        """
        bowtie2-build --threads {threads} {input.reference} {params.bowtie2_index} > {log} 2>&1
        bowtie2 {BOWTIE2_SE_OPTIONS} {BOWTIE2_OPTIONS} --threads {threads} -x {params.bowtie2_index} -U {input.reads} -S {output} >> {log} 2>&1
        rm {params.bowtie2_index}*
        """

rule samtools_view:
    input:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}" + SAMTOOLS_VIEW_INPUT_SUFFIX
    output:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}" + SAMTOOLS_VIEW_OUTPUT_SUFFIX
    log:
        "logs/samtools-view-" + sample_info.get_reference_prefix() + "-{sample}-{technology}.log"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools view --threads {threads} -bS {input} -o {output} > {log} 2>&1
        rm {input}
        """

rule samtools_sort:
    input:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}" + SAMTOOLS_SORT_INPUT_SUFFIX
    output:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}" + SAMTOOLS_SORT_OUTPUT_SUFFIX
    log:
        "logs/samtools-sort-" + sample_info.get_reference_prefix() + "-{sample}-{technology}.log"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools sort --threads {threads} {input} -o {output} > {log} 2>&1
        rm {input}
        """

# This operates on sorted BAM files
rule samtools_merge_nanopore:
    input:
        lambda wildcards: [ "aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + wildcards.sample + "/nanopore_reads." + s + SAMTOOLS_MERGE_INPUT_SUFFIX for s in sample_info.get_nanopore_read_base_list(wildcards.sample)]
    output:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/nanopore_reads_merged" + SAMTOOLS_MERGE_OUTPUT_SUFFIX
    log:
        "logs/samtools-merge-" + sample_info.get_reference_prefix() + "-{sample}-nanopore.log"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools merge --threads {threads} {input} -o {output} > {log} 2>&1
        rm {input}
        """
    
rule samtools_merge_illumina:
    input:
        lambda wildcards: [ "aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + wildcards.sample + "/illumina_reads." + s + ".SE" + SAMTOOLS_MERGE_INPUT_SUFFIX for s in sample_info.get_illumina_SE_read_base_list(wildcards.sample)] + [ "aligned-reads-" + sample_info.get_reference_prefix() + "/data/" + wildcards.sample + "/illumina_reads." + s + ".PE" + SAMTOOLS_MERGE_INPUT_SUFFIX for s in sample_info.get_illumina_PE_read_base_list(wildcards.sample)]
    output:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/illumina_reads_merged" + SAMTOOLS_MERGE_OUTPUT_SUFFIX
    log:
        "logs/samtools-merge-" + sample_info.get_reference_prefix() + "-{sample}-illumina.log"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools merge --threads {threads} {input} -o {output} > {log} 2>&1
        rm {input}
        """

rule samtools_bam_index:
    input:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam"
    output:
        "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/{technology}.bam.bai"
    log:
        "logs/samtools-index-" + sample_info.get_reference_prefix() + "-{sample}-{technology}.log"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools index -@ {threads} {input} > {log} 2>&1
        """

rule copy_convert_references:
    input:
        references = lambda wildcards: sample_info.get_reference_list(wildcards.sample)
    output:
        fasta = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta",
        gff3 = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.gff3",
    log:
        "logs/copy-convert-references-" + sample_info.get_reference_prefix() + "-{sample}.log"
    params:
        reference_arguments = lambda wildcards: sample_info.get_reference_arguments(wildcards.sample)
    conda:
        "../envs/breseq.yml"
    shell:
        """
        breseq CONVERT-REFERENCE -o {output.fasta} -f FASTA {params.reference_arguments} > {log} 2>&1
        breseq CONVERT-REFERENCE -o {output.gff3} --no-sequence -f GFF3 {params.reference_arguments} >> {log} 2>&1
        """

rule samtools_fasta_index:
    input:
        fasta = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta"
    output:
        fai = "aligned-reads-" + sample_info.get_reference_prefix() + "/data/{sample}/reference.fasta.fai"
    log:
        "logs/samtools-faidx-" + sample_info.get_reference_prefix() + "-{sample}.log"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools faidx {input.fasta} > {log} 2>&1
        """
