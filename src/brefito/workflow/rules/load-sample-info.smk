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

class SampleInfo():
    file_info = []

    sample_set = set({})
    sample_list = []

    # Dictionaries with keys that are samples and entries that are download (local) paths
    reference_lists = {}
    illumina_SE_read_lists = {}
    illumina_R1_read_lists = {}
    illumina_R2_read_lists = {}
    nanopore_read_lists = {}

    # used to check for duplicates and deconflict
    used_local_paths_set = set({})
    used_remote_paths_set = set({})

    # The class "constructor" - It's actually an initializer 
    def __init__(self, sample_info_csv_name):

        with open(sample_info_csv_name, encoding='utf-8-sig') as data_file:
            data_reader = csv.DictReader(data_file, delimiter=',', quotechar='"')
            for row in data_reader:
                print(row)
                file_name = os.path.basename(row['path'])
                base_file_name = os.path.basename(row['path'])
                simplified_read_name = self.get_simplified_read_name(base_file_name)
                file_name_without_extension = os.path.splitext(base_file_name)[0]

                if (row['type'] == "nanopore"):
                    local_path = os.path.join("nanopore_reads", simplified_read_name + ".fastq.gz")
                    self.append({ 
                        'sample' : row['sample'],
                        'type' : "nanopore", 
                        'file_type' : "nanopore", 
                        'remote_path' : row['path'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina"):
                    local_path = os.path.join("illumina_reads", simplified_read_name + ".fastq.gz")
                    self.append({ 
                        'sample' : row['sample'],
                        'type' : "illumina-SE", 
                        'file_type' : "illumina", 
                        'remote_path' : row['path'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R1"):
                    local_path = os.path.join("illumina_reads", simplified_read_name + ".R1.fastq.gz")
                    self.append({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R1", 
                        'file_type' : "illumina", 
                        'remote_path' : row['path'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R2"):
                    local_path = os.path.join("illumina_reads", simplified_read_name + ".R2.fastq.gz")
                    self.append({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R2", 
                        'file_type' : "illumina", 
                        'remote_path' : row['path'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "reference"):
                    local_path = os.path.join("references", file_name)
                    self.append({ 
                        'sample' : row['sample'],
                        'type' : "reference", 
                        'file_type' : "reference", 
                        'remote_path' : row['path'],
                        'local_path' : local_path
                        })
                else:
                    print("Skipping row with unknown type:" + row['type'])

        self.sample_list = list(self.sample_set)

    ## Check to see if the user is going to clobber their own files
    ## by having multiple remote paths mapping to the same download path

    def append(self, row):

        # Do checks that lead to not adding to the list first
        if self.remote_path_is_duplicate(row['remote_path']):
            print("Duplicate remote path encountered. Skipping entry:")
            print(row)
            return

        self.sample_set.add(row['sample'])
        row['local_path'] = self.deconflict_local_path(row['local_path'])
        self.file_info.append(row)

        print(row)

        sample = row['sample']
        if row['type'] == "reference":
            if not sample in self.reference_lists.keys():
                self.reference_lists[sample] = []
            self.reference_lists[sample].append(row['local_path'])
        elif row['type'] == "nanopore":
            if not sample in self.nanopore_read_lists.keys():
                self.nanopore_read_lists[sample] = []
            self.nanopore_read_lists[sample].append(row['local_path'])
        elif row['type'] == "illumina":
            if not sample in self.illumina_read_lists.keys():
                self.illumina_SE_read_lists[sample] = []
            self.illumina_SE_read_lists[sample].append(row['local_path'])
        elif row['type'] == "illumina-R1":
            if not sample in self.illumina_read_lists.keys():
                self.illumina_R1_read_lists[sample] = []
            self.illumina_R1_read_lists[sample].append(row['local_path'])
        elif row['type'] == "illumina-R2":
            if not sample in self.illumina_read_lists.keys():
                self.illumina_R2_read_lists[sample] = []
            self.illumina_R2_read_lists[sample].append(row['local_path'])

    def deconflict_local_path(self, this_local_path):
        new_local_path = this_local_path
        if (new_local_path in self.used_local_paths_set):
            i=1
            while new_local_path + "_" + str(i) in self.used_local_paths_set:
                i = i + 1
            new_local_path = new_local_path + "_" + str(i)
        self.used_local_paths_set.add(new_local_path)
        return new_local_path


    def remote_path_is_duplicate(self, this_remote_path):
        if (this_remote_path in self.used_remote_paths_set):
            return True
        self.used_remote_paths_set.add(this_remote_path)
        return False

    ## We want the read names to be standardized... this should do it in most cases
    def get_simplified_read_name(self, in_read_name):
        new_read_name = in_read_name
        new_read_name = new_read_name.replace("_R1", "").replace(".R1", "")
        new_read_name = new_read_name.replace("_R2", "").replace(".R2", "")
        new_read_name = new_read_name.replace(".fastq", "").replace(".gz", "")
        return new_read_name

    def get_nanopore_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.nanopore_read_lists[sample]])

    def get_nanopore_read_list(self, sample):
        return self.nanopore_read_lists[sample]

    def get_illumina_SE_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.illumina_SE_read_lists[sample]])
    
    def get_illumina_PE_read_arguments(self, sample, argument_R1_prefix='', argument_R2_prefix=''):
        args = ""
        if length(self.illumina_R1_read_lists[sample]) != length(self.illumina_R2_read_lists[sample]):
            print("Sample does not have an equal number of R1 and R2 entries: " + sample)
            sys.exit(1)

        arg_list = []
        for i in range(1,length(self.illumina_R1_read_lists[sample])):
            arg_list.append(argument_R1_prefix + self.illumina_R1_read_lists[sample][i] + argument_R2_prefix + self.illumina_R2_read_lists[sample][i])

        return " ".join([argument_prefix + item for item in self.reference_lists[sample]])

    def get_reference_arguments(self, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.reference_lists[sample]])

    def get_local_path_list(self):
        return [ d['local_path'] for d in self.file_info ]

    def get_remote_path_from_local_path(self, this_local_path):
        print(this_local_path)
        return [d['remote_path'] for d in self.file_info if d.get('local_path') == this_local_path]

    def get_sample_list(self):
        return self.sample_list

#### Initialize from config values
print(config)

data_csv_name = 'data.csv'
if "data_csv" in config:
    data_csv_name = config["data_csv"]

print("Loading Sample Info from file: " + data_csv_name)
sample_info = SampleInfo(data_csv_name)

print("CONFIG VARIABLES")

REMOTE_BASE_PATH='.'
if "remote_path" in config:
    REMOTE_BASE_PATH = config["remote_path"]
print("  remote_path=" + REMOTE_BASE_PATH)

BOOKMARK = ""
if "bookmark" in config.keys():
    BOOKMARK = config["bookmark"]
print("  bookmark=" + BOOKMARK)




