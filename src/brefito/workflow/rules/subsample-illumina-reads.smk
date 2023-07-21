include: "trim-illumina-reads.smk"

READ_NUMS = ["1", "2"]

rule subsample_illumina_reads:
    input:
        expand("01_trimmed_illumina_reads/{{dataset}}.R{read_num}.fastq.gz", read_num=READ_NUMS)
    output:
        expand("03_subsampled_illumina_reads/{{dataset}}/sample_{{assembly_id}}_R{read_num}.fastq.gz", read_num=READ_NUMS)
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        NOMINALDEPTH=80
        seqkit stats {input} -T > 01_trimmed_reads/{wildcards.dataset}.short_reads.stats.tsv
        BASES=`tail -n -2 01_trimmed_reads/{wildcards.dataset}.short_reads.stats.tsv | echo $(( $( cut -d$'\t' -f5 | paste -s -d+ - ) ))`
        echo $BASES
        READS=`tail -n -2 01_trimmed_reads/{wildcards.dataset}.short_reads.stats.tsv | echo $(( $( cut -d$'\t' -f4 | paste -s -d+ - ) ))`
        echo $READS
        READPROPORTION=`python -c "print(min(1, $BASES / {config[genome_size]} * $NOMINALDEPTH))"`
        echo $READPROPORTION

        MYRANDSEED=$RANDOM
        zcat {input[0]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output[0]}
        zcat {input[1]} | seqkit sample -s $MYRANDSEED -p $READPROPORTION -o {output[1]}
        """