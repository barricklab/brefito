## NOTES
# Spades errors out on newer mac operating systems and new spades versions with fixes
# can't be installed through conda. (Ugh)

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "trim-illumina-reads.smk"
include: "filter-nanopore-reads.smk"

rule subsample_nanopore_reads:
    input:
        "nanopore-reads0filtered/{sample}.fastq"
    output:
        "nanopore-reads-subsampled/{sample}.fastq"
    conda:
        "../envs/trycycler.yml"
    threads: 1
    shell:
        """
        trycycler subsample --genome_size {config[genome_size]} --reads {input} --count 2 --out_dir nanopore-reads-subsampled/{wildcards.sample}
        mv nanopore-reads-subsampled/{wildcards.sample}/sample_01.fastq {output}
        rm -r nanopore-reads-subsampled/{wildcards.sample}
        """

rule gzip_subsampled_nanopore_reads:
    input:
        "nanopore-reads-subsampled/{sample}.fastq"
    output:
        "nanopore-reads-subsampled/{sample}.fastq.gz"
    threads: 1
    shell:
        """
        gzip {input}
        """

READ_NUMS = ["1", "2"]

rule subsample_illumina_reads:
    input:
        expand("illumina-reads-trimmed/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        reads = expand("illumina-reads-subsampled/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    log:
        "logs/subsample-illumina-reads/{dataset}.log"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        echo "subsample-illumina-reads" > {log}
        echo "Provided assembly size: {config[genome_size]}" >> {log}
        NOMINALDEPTH=80
        echo "Target nominal read depth: $NOMINALDEPTH\n" >> {log}
        seqkit stats {input} -T > illumina-reads-subsampled/{wildcards.dataset}.short_reads.stats.tsv
        BASES=`tail -n -2 illumina-reads-subsampled/{wildcards.dataset}.short_reads.stats.tsv | awk '{{sum += $5}} END {{print sum}}'`
        echo "Total bases: $BASES" >> {log}
        READS=`tail -n -2 illumina-reads-subsampled/{wildcards.dataset}.short_reads.stats.tsv | awk '{{sum += $4}} END {{print sum}}'`
        echo "Total reads: $READS" >> {log}
        READPROPORTION=`echo "scale=9; $BASES / {config[genome_size]} * $NOMINALDEPTH" | bc`
        READPROPORTION=$(awk -v val="$READPROPORTION" 'BEGIN {{if (val > 1) print 1; else print val}}')
        echo "Calculated read proportion to sample: $READPROPORTION" >> {log}

        MYRANDSEED=$RANDOM
        echo "Random seed: $MYRANDSEED" >> {log}

        gunzip -c {input[0]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output.reads[0]} 2>> {log}
        gunzip -c {input[1]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output.reads[1]} 2>> {log}
        """

# Helper functions for setting up arguments
def find_available_read_files(wildcards):
    nanopore_files = [os.path.join("nanopore-reads-subsampled", os.path.basename(d)) for d in sample_info.get_nanopore_read_list(wildcards.sample)]
    illumina_files = [os.path.join("illumina-reads-subsampled", os.path.basename(d)) for d in sample_info.get_illumina_read_list(wildcards.sample)]
    return illumina_files + nanopore_files


def get_paired_short_read_args(wildcards):
    return sample_info.get_illumina_PE_read_arguments(wildcards.sample, "-1 illumina-reads-subsampled/", "-2 illumina-reads-subsampled/")

def get_long_read_args(wildcards):
    return sample_info.get_nanopore_read_arguments(wildcards.sample, "-l ")

rule assemble_all:
    input:
        ["assemblies/" + s + ".fasta" for s in sample_info.get_sample_list()]
    default_target: True

rule assemble_with_unicycler:
    input:
        reads = lambda wildcards: find_available_read_files(wildcards)
    output:
        fasta = "assemblies/{sample}.fasta",
        gfa = "assemblies/{sample}.gfa",
        directory = temp(directory("unicycler-assembly/{sample}"))
    log:
        "logs/assemble_with_unicycler/{sample}.log"
    params:
        paired_short_read_args =  lambda wildcards: get_paired_short_read_args(wildcards),
        long_read_args = lambda wildcards: get_long_read_args(wildcards)
    conda:
        "../envs/unicycler.yml"
    threads: 3
    shell:
        """
        mkdir -p assemblies
        unicycler --threads {threads} {params.paired_short_read_args} {params.long_read_args} -o {output.directory} > {log} 2>&1
        mv {output.directory}/assembly.fasta {output.fasta}
        mv {output.directory}/assembly.gfa {output.gfa}
        """