# CONFIG OPTIONS
#
# bookmark = bookmark (REQUIRED)
# data_csv = path to CSV file descring(default = 'data.csv')
# remote_path = (default = '.')
#
# DATA_CSV
#
# This file must have the following columns:
#  1 sample = main sample name, must match for reads that will be used together
#  2 type = either 'nanopore', 'illumina_R1', 'illumina_R2'
#  3 path = path for file on remote relative to remote_path
#
# SETUP BOOKMARK
#
# In order to connect you need to set up an lftp bookmark. 
# One way to do this is to run a command like this within your lftp shell
# bookmark utbox ftps://'<username>':<password>@ftp.box.com:990
# SECURITY WARNING: this creates a file - to which only you have access
# by default - which has your password as plain text under ~/.local
# if you provide a password

import csv

print(config)
BOOKMARK = ""
if "bookmark" in config.keys():
    BOOKMARK = config["bookmark"]

file_info = []
URL_dict = {}

data_csv = 'data.csv'
if "data_csv" in config:
    data_csv = config["data_csv"]

REMOTE_BASE_PATH='.'
if "remote_path" in config:
    REMOTE_BASE_PATH = config["remote_path"]

with open(data_csv) as data_file:
    data_reader = csv.DictReader(data_file, delimiter=',', quotechar='"')
    for row in data_reader:
        if (row['type'] == "nanopore"):
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "nanopore", 
                'remote_path' : row['path'],
                'download_path' : "nanopore_reads/" + row['sample'] + ".fastq.gz"
                })
            URL_dict["nanopore-" + row['sample']] = row['path']
        elif (row['type'] == "illumina_R1"):
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "illumina", 
                'remote_path' : row['path'],
                'download_path' : "illumina_reads/" + row['sample'] + ".R1.fastq.gz"
                })
            URL_dict["illumina-" + row['sample'] + ".R1"] = row['path']
        elif (row['type'] == "illumina_R2"):
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "illumina", 
                'remote_path' : row['path'],
                'download_path' : "illumina_reads/" + row['sample'] + ".R2.fastq.gz"
                })
            URL_dict["illumina-" + row['sample'] + ".R2"] = row['path']
        else:
            print("Skipping row with unknown type:" + row['type'])

print(file_info)
print(URL_dict)

def download_all():
    print([ d['download_path'] for d in file_info])
    return([ d['download_path'] for d in file_info])

rule all:
    input: download_all()

rule download_reads:
    output:
        output_path = "{read_type}_reads/{sample}.fastq.gz",
    params:
        URL=lambda wildcards: URL_dict[wildcards.read_type + "-" + wildcards.sample],
        bookmark = BOOKMARK,
        download_path = "{read_type}_reads/temp_{sample}.fastq.gz",
        lftp_commands_file = "{read_type}_reads/temp_{sample}.lftp.commands"
    threads: 1
    resources:
        # This is an invented resource to prevent opening too many download connections at once!
        connections=1
    shell:
        """
        echo 'cd "{REMOTE_BASE_PATH}"' > {params.lftp_commands_file}
        echo 'get "{params.URL}" -o "{params.download_path}"' >> {params.lftp_commands_file} 
        #cat {params.lftp_commands_file} 
        lftp {params.bookmark} < {params.lftp_commands_file}
        rm {params.lftp_commands_file} 
        mv {params.download_path} {output.output_path}
        """    
