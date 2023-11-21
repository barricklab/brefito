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
import os

print(config)
BOOKMARK = ""
if "bookmark" in config.keys():
    BOOKMARK = config["bookmark"]

file_info = []
URL_dict = {}

reference_lists = {}
illumina_read_lists = {}
nanopore_read_lists = {}
sample_set = set({})
sample_list = []

data_csv = 'data.csv'
if "data_csv" in config:
    data_csv = config["data_csv"]

REMOTE_BASE_PATH='.'
if "remote_path" in config:
    REMOTE_BASE_PATH = config["remote_path"]

## Check to see if the user is going to clobber their own files
## by having multiple remote paths mapping to the same download path
used_download_paths = {}
def check_for_duplicate_download_paths(download_path, remote_path):
    if download_path in used_download_paths.keys():
        if (used_download_paths[download_path] != remote_path):
            print("ERROR: Multiple remote paths map to the same download path!")
            print("  " + used_download_paths[download_path] )
            print("  " + remote_path)
            sys.exit(1)
    else:
        used_download_paths[download_path] = remote_path

## We want the read names to be standardizes... this should do it in most cases
def get_simplified_read_name(in_read_name):
    new_read_name = in_read_name
    new_read_name = new_read_name.replace("_R1", "").replace(".R1", "")
    new_read_name = new_read_name.replace("_R2", "").replace(".R2", "")
    new_read_name = new_read_name.replace(".fastq", "").replace(".gz", "")
    return new_read_name

with open(data_csv, encoding='utf-8-sig') as data_file:
    data_reader = csv.DictReader(data_file, delimiter=',', quotechar='"')
    for row in data_reader:
        #print(row)
        file_name = os.path.basename(row['path'])
        base_file_name = os.path.basename(row['path'])
        simplified_read_name = get_simplified_read_name(base_file_name)
        file_name_without_extension = os.path.splitext(base_file_name)[0]

        sample_set.add(row['sample'])

        if (row['type'] == "nanopore"):
            download_path = os.path.join("nanopore_reads", simplified_read_name + ".fastq.gz")
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "nanopore", 
                'remote_path' : row['path'],
                'download_path' : download_path
                })
            URL_dict[download_path] = row['path']
            check_for_duplicate_download_paths(download_path, row['path'])

        elif (row['type'] == "illumina_R1"):
            download_path = os.path.join("illumina_reads", simplified_read_name + ".R1.fastq.gz")
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "illumina", 
                'remote_path' : row['path'],
                'download_path' : download_path
                })
            URL_dict[download_path] = row['path']
            check_for_duplicate_download_paths(download_path, row['path'])


        elif (row['type'] == "illumina_R2"):
            download_path = os.path.join("illumina_reads", simplified_read_name + ".R2.fastq.gz")
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "illumina", 
                'remote_path' : row['path'],
                'download_path' : download_path
                })
            URL_dict[download_path] = row['path']
            check_for_duplicate_download_paths(download_path, row['path'])

        elif (row['type'] == "reference"):
            download_path = os.path.join("references", file_name)
            file_info.append({ 
                'sample' : row['sample'],
                'file_type' : "reference", 
                'remote_path' : row['path'],
                'download_path' : download_path
                })
            URL_dict[download_path] = row['path']
            check_for_duplicate_download_paths(download_path, row['path'])
        else:
            print("Skipping row with unknown type:" + row['type'])

    sample_list = list(sample_set)
    for line in file_info:
        sample = line['sample']
        if line['file_type'] == "reference":
            if not sample in reference_lists.keys():
                reference_lists[sample] = []
            reference_lists[sample].append(line['download_path'])
        elif line['file_type'] == "nanopore":
            if not sample in nanopore_read_lists.keys():
                nanopore_read_lists[sample] = []
            nanopore_read_lists[sample].append(line['download_path'])
        elif line['file_type'] == "illumina":
            if not sample in illumina_read_lists.keys():
                illumina_read_lists[sample] = []
            illumina_read_lists[sample].append(line['download_path'])