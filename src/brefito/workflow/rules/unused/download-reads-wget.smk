<NONFUNCTIONAL>

# Config options
#
# data_csv = path to CSV file descring(default = "data.csv")
#
# The data_csv file must have the following columns:
#  1 sample
#  2 type = either nanopore or illumina
#  3 path = path on remove relative to base
# sampletype of reads (nanopore or illumina), path on remote relative to base


# Must set --config url="https://something.something/path/to/base/"


#
"data.csv"

import csv

#dictionary
nanopore_URL_lookup = dict()
illumina_R1_URL_lookup = dict()
illumina_R2_URL_lookup = dict()

with open('data.csv') as data_file:

    data_reader = csv.DictReader(data_file, delimiter=',', quotechar='"')

    for row in data_reader:
        if (row['type'].lower() == "nanopore"):
            nanopore_URL_lookup[row['sample']] = row['path']
        elif (row['type'].lower() == "illumina_R1"):
            illumina_R1_URL_lookup[row['sample']] = row['path']
        elif (row['type'].lower() == "illumina_R2"):
            illumina_R2_URL_lookup[row['sample']] = row['path']
        else:
            print("Skipping row with unknown type:" + row['type'])

print(type(nanopore_URL_lookup))
print(nanopore_URL_lookup)

def download_all():
    nanopore_read_files = [ "nanopore_reads/" + d + ".fastq.gz" for d in nanopore_URL_lookup.keys()]
    illumina_read_R1_files = [ "illumina_reads/" + d + ".R1.fastq.gz" for d in illumina_R1_URL_lookup.keys()]
    illumina_read_R2_files = [ "illumina_reads/" + d + ".R2fastq.gz" for d in illumina_R2_URL_lookup.keys()]

    print(nanopore_read_files)
    return (nanopore_read_files + illumina_read_R1_files + illumina_read_R2_files)

def find_URL(wildcards):
    clusters = [str(x) for x in glob.glob(f"05_trycycler/{wildcards.dataset}/cluster_*/")]
    clusters = list(map(lambda x: x + "8_medaka.fasta", clusters))
    print(clusters)
    return sorted(clusters) 

rule all:
    input: download_all()

rule download_nanopore_reads:
    output:
        output_path = "nanopore_reads/{sample}.fastq.gz",
        download_path = temp("downloads/nanopore_reads/{sample}.fastq.gz")
    log: 
        "logs/evaluate_breseq_nanopore_reads_{sample}.log"
    params:
        URL=lambda wildcards: nanopore_URL_lookup[wildcards.sample]
    threads: 1
    shell:
        """
        wget -O {output.download_path} {config[url]}/{params.URL}
        mv {output.download_path} {output.output_path}
        """    

rule download_illumina_R1_reads:
    output:
        output_path = "nanopore_reads/{sample}.R1.fastq.gz",
        download_path = temp("downloads/nanopore_reads/{sample}.R1.fastq.gz")
    log: 
        "logs/evaluate_breseq_nanopore_reads_{sample}.log"
    params:
        URL=lambda wildcards: illumina_R1_URL_lookup[wildcards.sample]
    threads: 1
    shell:
        """
        wget -O {output.download_path} {config[url]}/{params.URL}
        mv {output.download_path} {output.output_path}
        """    

rule download_illumina_R2_reads:
    output:
        output_path = "nanopore_reads/{sample}.R1.fastq.gz",
        download_path = temp("downloads/nanopore_reads/{sample}.R2.fastq.gz")
    log: 
        "logs/evaluate_breseq_nanopore_reads_{sample}.log"
    params:
        URL=lambda wildcards: illumina_R2_URL_lookup[wildcards.sample]
    threads: 1
    shell:
        """
        wget -O {output.download_path} {config[url]}/{params.URL}
        mv {output.download_path} {output.output_path}
        """    