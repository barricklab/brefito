try: sample_info
except NameError: 
    include: "load-sample-info.smk"

include: "filter-nanopore-reads.smk"

rule subsample_nanopore_reads:
    input:
        "nanopore_reads_filtered/{sample}.fastq"
    output:
        "nanopore_reads_subsampled/{sample}.fastq"
    conda:
        "../envs/trycycler.yml"
    threads: 1
    shell:
        """
        trycycler subsample --genome_size {config[genome_size]} --reads {input} --count 2 --out_dir nanopore_reads_subsampled/{wildcards.sample}
        mv nanopore_reads_subsampled/{wildcards.sample}/sample_01.fastq {output}
        rm -r nanopore_reads_subsampled/{wildcards.sample}
        """

rule gzip_subsampled_nanopore_reads:
    input:
        "nanopore_reads_subsampled/{sample}.fastq"
    output:
        "nanopore_reads_subsampled/{sample}.fastq.gz"
    threads: 1
    shell:
        """
        gzip {input}
        """

READ_NUMS = ["1", "2"]

rule subsample_illumina_reads:
    input:
        expand("illumina_reads_trimmed/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("illumina_reads_subsampled/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        NOMINALDEPTH=80
        seqkit stats {input} -T > illumina_reads_subsampled/{wildcards.dataset}.short_reads.stats.tsv
        BASES=`tail -n -2 illumina_reads_subsampled/{wildcards.dataset}.short_reads.stats.tsv | awk '{{sum += $5}} END {{print sum}}'`
        echo $BASES
        READS=`tail -n -2 illumina_reads_subsampled/{wildcards.dataset}.short_reads.stats.tsv | awk '{{sum += $4}} END {{print sum}}'`
        echo $READS
        READPROPORTION=`echo "scale=9; $BASES / {config[genome_size]} * $NOMINALDEPTH" | bc`
        READPROPORTION=$(awk -v val="$READPROPORTION" 'BEGIN {{if (val > 1) print 1; else print val}}')
        echo $READPROPORTION

        MYRANDSEED=$RANDOM
        gzcat {input[0]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output[0]}
        gzcat {input[1]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output[1]}
        """

# Helper functions for setting up arguments
def find_available_read_files(wildcards):

    nanopore_files = [os.path.join("nanopore_reads_subsampled", os.path.basename(d)) for d in sample_info.get_nanopore_read_list(wildcards.sample)]
    illumina_files = [os.path.join("illumina_reads_subsampled", os.path.basename(d)) for d in sample_info.get_illumina_read_list(wildcards.sample)]

    print(illumina_files + nanopore_files)

    return illumina_files + nanopore_files


def get_paired_short_read_args(wildcards):
    print("illumina_read_args")
    args = sample_info.get_illumina_PE_read_arguments(wildcards.sample, "-1 illumina_reads_subsampled/", "-2 illumina_reads_subsampled/")
    print(args)
    return args

def get_long_read_args(wildcards):
    print("nanopore_read_args")
    args = sample_info.get_nanopore_read_arguments(wildcards.sample, "-l ")
    print(args)
    return args

rule assemble_all:
    input:
        ["assemblies/" + s + ".fasta" for s in sample_info.get_sample_list()]
    default_target: True

rule assemble_with_unicycler:
    priority: 1
    input:
        reads = lambda wildcards: find_available_read_files(wildcards)
    output:
        fasta = "assemblies/{sample}.fasta",
        gfa = "assemblies/{sample}.gfa",
        directory = temp(directory("unicycler_assembly/{sample}"))
    params:
        paired_short_read_args =  lambda wildcards: get_paired_short_read_args(wildcards),
        long_read_args = lambda wildcards: get_long_read_args(wildcards)
    conda:
        "../envs/unicycler.yml"
    threads: 24
    shell:
        """
        mkdir -p {output.directory}
        unicycler --threads {threads} {params.paired_short_read_args} {params.long_read_args} -o {output.directory}
        cp {output.directory}/assembly.fasta {output.fasta}
        cp {output.directory}/assembly.gfa {output.gfa}
        rm -r {output.directory}
        """