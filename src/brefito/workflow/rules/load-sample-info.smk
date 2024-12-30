# CONFIG OPTIONS
#
# bookmark = bookmark (REQUIRED)
# data_csv = path to CSV file (default = 'data.csv')
# remote_path = (default = '.')
#
# DATA_CSV
#
# This file must have the following columns:
#  1 sample = main sample name, must match for reads that will be used together
#  2 type = either a file type: 'nanopore', 'illumina', 'illumina_R1', 'illumina_R2'
#           or an option setting, like breseq_option
#  3 setting = for file types, a way of downloading the file:
#                 lftp@bookmark://path/to/file (if bookmark is present and refers to a location)
#                 wget://example.com/path/to/file
#              for others, the value of the setting
#                 XXXXXXX
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
import glob
import copy

def rreplace(s, old, new, maxreplace):
    return(new.join(s.rsplit(old, maxreplace)))

def rreplace1(s, old, new):
    return(rreplace(s,old,new,1))

class SampleInfo():
    file_info = []

    sample_set = set({})
    sample_list = []

    # Dictionary of lists with keys that are types of data and then lists of files
    file_lists_by_sample_by_type = {}
    option_lists_by_sample_by_type = {}

    # used to check for duplicates and deconflict
    local_to_remote_path_mapping = {}

    # keeps track of all samples that use a local path
    local_paths_by_sample = {}

    reference_prefix = "references"

    # The class constructor
    def __init__(self, sample_info_csv_name=None):
        if sample_info_csv_name != None:
            self.init_from_CSV_file(sample_info_csv_name)
        else:
            self.init_from_directory_structure()

    def init_from_CSV_file(self, sample_info_csv_name=None):

        with open(sample_info_csv_name, encoding='utf-8-sig') as data_file:
            data_reader = csv.DictReader(data_file, delimiter=',', quotechar='"')
            for row in data_reader:
                #print(row)

                #Ignore commented rows, where first item begins with #
                if 'sample' in row:
                    if row['sample'][0]=='#':
                        continue

                # Add 'path' for 'setting' for backwards compatibility
                if 'path' in row:
                    if ('setting' in row):
                        print("Both 'path' and 'setting' found for row. Only using 'setting'.")
                    else:
                        row["setting"] = row["path"]

                # Copy over certain paths to allow for some choices.
                file_name = os.path.basename(row['setting'])
                base_file_name = os.path.basename(row['setting'])

                #And then also remove .fastq here
                simplified_read_file_base_name = self.get_simplified_read_file_base_name(base_file_name)

                if (row['type'] == "nanopore"):
                    local_path = os.path.join("nanopore-reads", simplified_read_file_base_name + ".fastq.gz")
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "nanopore", 
                        'file_type' : "nanopore", 
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina") or (row['type'] == "illumina-SE"):
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".SE.fastq.gz")
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "illumina-SE", 
                        'file_type' : "illumina", 
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R1"):
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".R1.fastq.gz")
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R1", 
                        'file_type' : "illumina", 
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R2"):
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".R2.fastq.gz")
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R2", 
                        'file_type' : "illumina", 
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-pair") or (row['type'] == "illumina-paired") or (row['type'] == "illumina-PE"):
                    #Must have {1/2} in the name
                    if not "{1|2}" in row['setting']:
                        print("Did not find {1|2} in name of illumina-PE entry:" + row['setting'])
                        print("Skipping this entry.")
                        continue

                    # Don't do the automatic read name simplification step
                    simplified_read_file_base_name_1 = self.get_simplified_read_file_base_name(base_file_name.replace("{1|2}", "1"))
                    simplified_read_file_base_name_2 = self.get_simplified_read_file_base_name(base_file_name.replace("{1|2}", "2"))

                    remote_path_1 = row['setting'].replace("{1|2}", "1")
                    remote_path_2 = row['setting'].replace("{1|2}", "2")

                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R1",
                        'file_type' : "illumina", 
                        'remote_path' : remote_path_1,
                        'local_path' : os.path.join("illumina-reads", simplified_read_file_base_name_1 + ".R1.fastq.gz")
                        })
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "illumina-R2", 
                        'file_type' : "illumina", 
                        'remote_path' : remote_path_2,
                        'local_path' : os.path.join("illumina-reads", simplified_read_file_base_name_2 + ".R2.fastq.gz")
                        })

                elif (row['type'] == "reference"):
                    local_path = os.path.join("references", file_name)
                    self.add_file({ 
                        'sample' : row['sample'],
                        'type' : "reference", 
                        'file_type' : "reference", 
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })
                else:
                    # add this to the options dictionary
                    if not row['sample'] in self.option_lists_by_sample_by_type:
                        self.option_lists_by_sample_by_type[row['sample']] = {}
                    if not row['type'] in self.option_lists_by_sample_by_type[row['sample']]:
                        self.option_lists_by_sample_by_type[row['sample']][row['type']] = []
                    self.option_lists_by_sample_by_type[row['sample']][row['type']].append(row['setting'])

                #else:
                #    print("Skipping row with unknown type:" + row['type'])

        self.sample_list = list(self.sample_set)
        return True

    # Helper function for identifying valid input files
    # Returns a dict with key as sample name determined from file
    def find_matching_input_files(self, in_base_path, *in_file_endings):

        matching_input_files = {}

        for in_file_ending in in_file_endings:

            existing_files = glob.glob(os.path.join(in_base_path, "*" + in_file_ending))

            #print(os.path.join(base_path,"input", "*."+file_ending))
            for this_input_file in existing_files:
                this_file_name=os.path.basename(this_input_file)

                #print(this_file_name)
                short_name = re.findall("^(.+)" + re.escape(in_file_ending) + "$", this_file_name)
                matching_input_files[short_name[0]] = this_input_file

        return(matching_input_files)
        

    def init_from_directory_structure(self):
         # What files are available?

        assemblies_path = "assemblies"
        references_path = "references"
        nanopore_input_path = "nanopore-reads"
        illumina_input_path = "illumina-reads"

        print()
        print("Nanopore read files found (*.fastq.gz) in " + nanopore_input_path)
        print()

        input_nanopore_files = self.find_matching_input_files(nanopore_input_path, ".fastq.gz")
        for k in input_nanopore_files: print("    " + k)
        if (len(input_nanopore_files.items()) == 0) : print("    " + "NONE FOUND")
        for k, i in input_nanopore_files.items(): 
            self.add_file({ 
                        'sample' : k,
                        'type' : "nanopore", 
                        'file_type' : "nanopore", 
                        'remote_path' : None,
                        'local_path' : i
                        })

        print()
        print("Paired-end Illumina read files found (*R(1|2)*.fastq.gz) in " + illumina_input_path)
        print()

        input_illumina_1_files = self.find_matching_input_files(illumina_input_path, ".R1.fastq.gz")
        input_illumina_2_files = self.find_matching_input_files(illumina_input_path, ".R2.fastq.gz")

        for key in set( list(input_illumina_1_files.keys()) + list(input_illumina_2_files.keys()) ):
            if key in input_illumina_1_files.keys(): print ("    " + str(key) + " : " + input_illumina_1_files[key])
            if key in input_illumina_2_files.keys(): print ("    " + str(key) + " : " + input_illumina_2_files[key])

        input_illumina_PE_file_names = set({})
        for k, i in input_illumina_1_files.items(): 
            self.add_file({ 
                        'sample' : k,
                        'type' : "illumina-R1", 
                        'file_type' : "illumina", 
                        'remote_path' : None,
                        'local_path' : i
                        })
            input_illumina_PE_file_names.add(i)

        for k, i in input_illumina_2_files.items(): 
            self.add_file({ 
                        'sample' : k,
                        'type' : "illumina-R2", 
                        'file_type' : "illumina", 
                        'remote_path' : None,
                        'local_path' : i
                        })
            input_illumina_PE_file_names.add(i)



        if (len(input_illumina_1_files.items()) == 0) and (len(input_illumina_2_files.items()) == 0): 
            print("    " + "NONE FOUND")

        print()
        print("Unpaired Illumina read files found (*.fastq.gz) in " + illumina_input_path)
        print()

        input_illumina_SE_files = self.find_matching_input_files(illumina_input_path, ".fastq.gz")
        for k, i in input_illumina_SE_files.items():
            if i in input_illumina_PE_file_names:
                continue
            self.add_file({ 
                'sample' : k,
                'type' : "illumina-SE", 
                'file_type' : "illumina", 
                'remote_path' : None,
                'local_path' : i
                })   
            print("    " + k + " : " + i)



        print()
        print("Genome assembly files found (*.fasta) in " + assemblies_path)
        print()
        input_assembly_files = self.find_matching_input_files(assemblies_path, ".fasta", ".fna", ".fa")
        for (k, v) in input_assembly_files.items(): print("    " + k + " : " + v)
        if (len(input_assembly_files.items()) == 0) : print("    " + "NONE FOUND")

        for k, i in input_assembly_files.items(): 
            self.add_file({ 
                        'sample' : k,
                        'type' : "assembly", 
                        'file_type' : "fasta", 
                        'remote_path' : None,
                        'local_path' : i
                        })

        input_main_reference_assembly_files = glob.glob("reference.fasta")
        input_main_reference_assembly_file = None
        input_main_reference_assembly_file_status = "UNKNOWN";
        if (len(input_main_reference_assembly_files) == 1):
            input_main_reference_assembly_file = input_main_reference_assembly_files[0]
            input_main_reference_assembly_file_status = "FOUND"
        else:
            input_main_reference_assembly_file_status = "NOT FOUND"
        print()
        print("Main reference genome assembly file found (reference.fasta)? " + input_main_reference_assembly_file_status)
        print()

        print()
        print("Reference genome assembly files found (*.fasta, *.fna, *.fa, *.genbank, *.gbk, *.gb, *.gff, *.gff3) in " + references_path)
        print()

        input_reference_files = self.find_matching_input_files(references_path, ".fasta", ".fna", ".fa", ".genbank", ".gbk", ".gb", ".gff", ".gff3")
        for (k, v) in input_reference_files.items(): print("    " + k + " : " + v)
        if (len(input_reference_files.items()) == 0) : print("    " + "NONE FOUND")

        for k, i in input_reference_files.items(): 
            self.add_file({ 
                        'sample' : k,
                        'type' : "reference", 
                        'file_type' : "reference", 
                        'remote_path' : None,
                        'local_path' : i
                        })

    def add_file(self, row):

        # Do checks that lead to not adding to the list first
        #if row['remote_path']!= None and self.remote_path_is_duplicate(row['remote_path']):
        #    print("Duplicate remote path encountered. Skipping entry:")
        #    print(row)
        #    return

        # * means to apply to all previous entries
        if row['sample'] == '*':
            this_sample_list = list(self.sample_set)
        else:
            this_sample_list = [row['sample']]

        self.deconflict_paths(row)

        # This stores a full local path
        self.local_to_remote_path_mapping[row['local_path']] = row['remote_path']

        # This stores which local path belongs with which sample
        for this_sample in this_sample_list:
            if not this_sample in self.local_paths_by_sample.keys():
                self.local_paths_by_sample[this_sample] = []
            self.local_paths_by_sample[this_sample].append(row['local_path'])

        # Now we want a base path for the file_lists_by_sample_by_type dict
        row['local_path'] = os.path.split(row['local_path'])[1]

        for this_sample in this_sample_list:

            self.sample_set.add(this_sample)
            row['sample'] = this_sample

            #check for duplicates
            if row in self.file_info:
                print("Ignored duplicate sample info row:")
                print(row)
            else:
                self.file_info.append(copy.deepcopy(row))
                #print("Adding row")
                #print(row)

                if not this_sample in self.file_lists_by_sample_by_type.keys():
                    self.file_lists_by_sample_by_type[this_sample] = {}
                if not row['type'] in self.file_lists_by_sample_by_type[this_sample].keys():
                    self.file_lists_by_sample_by_type[this_sample][row['type']] = []
                self.file_lists_by_sample_by_type[this_sample][row['type']].append(row['local_path'])

    # Looks for matches in a list and splits them off (to allow multile extensions)
    def split_filename_and_extension(self, input_filename):
        current_file_name = input_filename
        extension_list = [""]
        valid_extension_set = set({".gz", ".fastq", ".fasta", ".fna", ".fa", ".genbank", ".gbk", ".gb", ".gff3", ".gff", ".SE", ".R1", ".R2"})
        valid_extension_found = True
        while valid_extension_found:
            current_file_name, next_file_extension = os.path.splitext(current_file_name)
            #print(current_file_name, "    ", next_file_extension)
            if next_file_extension in valid_extension_set:
                extension_list.append(next_file_extension)
            else:
                current_file_name = current_file_name + next_file_extension
                valid_extension_found = False

        return [current_file_name, "".join(reversed(extension_list)) ]

    #makes sure that different remote paths aren't mapped to the same local path
    def deconflict_paths(self, this_row):
        if this_row['local_path'] in self.local_to_remote_path_mapping.keys():
            if self.local_to_remote_path_mapping[this_row['local_path']] != this_row['remote_path']:
                i=1
                file_name, file_extension = self.split_filename_and_extension(this_row['local_path'])
                while file_name + "-" + str(i) + file_extension in self.local_to_remote_path_mapping.keys():
                    i = i + 1
                this_row['local_path'] = file_name + "-" + str(i) + file_extension

    ## We want the read names to be standardized... this should do it in most cases
    def get_simplified_read_file_base_name(self, in_read_name):
        
        no_ending_read_name = in_read_name

        # We remove the last extension no matter what it is. Usually it is *.gz
        # This may not be necessary or might even cause problem...
        no_ending_read_name = os.path.splitext(no_ending_read_name)[0]
        #new_read_name = rreplace1(new_read_name, ".gz", "")

        # Then, remove multiple line endings, but only of certain formats
        no_ending_read_name = rreplace1(no_ending_read_name, ".fastq", "")

        # Now, remove the first of these that we encounter. 
        # The idea is that if we do this to the file names for
        # each read in a pair, we will have the SAME name afterward

        # These conditionals make sure we only remove one instance!!
        new_read_name = no_ending_read_name

        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_R1.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_R2.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, ".R1.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, ".R2.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_R1_", "_")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_R2_", "_")
        if new_read_name==no_ending_read_name:
            if new_read_name.endswith("_R1"):
                new_read_name=new_read_name[:-3]
        if new_read_name==no_ending_read_name:
            if new_read_name.endswith("_R2"):
                new_read_name=new_read_name[:-3]
            if new_read_name.endswith(".R1"):
                new_read_name=new_read_name[:-3]
        if new_read_name==no_ending_read_name:
            if new_read_name.endswith(".R2"):
                new_read_name=new_read_name[:-3]
        if new_read_name==no_ending_read_name:
            if new_read_name.endswith("_1"):
                new_read_name=new_read_name[:-2]
        if new_read_name==no_ending_read_name:
            if new_read_name.endswith("_2"):
                new_read_name=new_read_name[:-2]
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_1.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_2.", ".")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_1_", "_")
        if new_read_name==no_ending_read_name:
            new_read_name = rreplace1(new_read_name, "_2_", "_")

        #print(new_read_name)
        return new_read_name

    def get_file_list(self, in_sample, in_type):
        #print(in_sample, in_type)
        if not in_sample in self.file_lists_by_sample_by_type.keys(): return []
        if not in_type in self.file_lists_by_sample_by_type[in_sample].keys(): return []
        return self.file_lists_by_sample_by_type[in_sample][in_type]

    def print_file_lists(self):

        print()
        main_reference_file_path = self.get_main_reference_file()
        if os.path.isfile(main_reference_file_path):
            print("MAIN REFERENCE FILE :: " + main_reference_file_path)
        else:
            print("MAIN REFERENCE FILE :: NOT FOUND")

        print()
        print("NUMBER OF SAMPLES :: " + str(len(self.get_sample_list())))

        for this_sample in sorted(self.get_sample_list()):
            print()
            print("SAMPLE :: " + this_sample)
            for this_type in sorted(self.file_lists_by_sample_by_type[this_sample].keys()):
                for this_item in self.file_lists_by_sample_by_type[this_sample][this_type]:
                    print("  " + this_type + " : " + this_item)
        #print(self.local_to_remote_path_mapping)

    def get_nanopore_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.get_file_list(sample, "nanopore")])

    def get_nanopore_read_list(self, sample):
        return self.get_file_list(sample, "nanopore")

    def get_nanopore_read_base_list(self, sample):
        return [os.path.split(n.replace(".fastq.gz", ""))[1] for n in self.get_file_list(sample, "nanopore")]

    def get_samples_with_nanopore_reads(self):
        return [sample for sample in self.get_sample_list() if len(self.get_nanopore_read_list(sample))]

    def get_illumina_read_list(self, sample):
        return self.get_file_list(sample, "illumina-SE") + self.get_file_list(sample, "illumina-R1") + self.get_file_list(sample, "illumina-R2")

    def get_illumina_SE_read_base_list(self, sample):
        return self.get_file_list(sample, "illumina-SE")

    def get_illumina_SE_read_list(self, sample):
        return [os.path.split(i_se.replace(".SE.fastq.gz", ""))[1] for i_se in self.get_file_list(sample, "illumina-SE")]

    def get_illumina_SE_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item] for item in self.get_file_list(sample, "illumina-SE"))
    
    def get_samples_with_illumina_SE_reads(self):
        return [sample for sample in self.get_sample_list() if len(self.get_illumina_SE_read_base_list(sample))]

    def get_illumina_PE_read_list(self, sample):
        return self.get_file_list(sample, "illumina-R1") + self.get_file_list(sample, "illumina-R2")

    def get_illumina_PE_read_base_list(self, sample):
        illumina_R1_read_list = sorted(self.get_file_list(sample, "illumina-R1"))
        illumina_R2_read_list = sorted(self.get_file_list(sample, "illumina-R2"))

        if len(illumina_R1_read_list) != len(illumina_R2_read_list):
            print("Sample does not have an equal number of R1 and R2 entries: " + sample)
            sys.exit(1)

        read_base_list = []
        for i in range(0,len(illumina_R1_read_list)):
            illumina_R1_base = illumina_R1_read_list[i].replace(".R1.fastq.gz", "")
            illumina_R2_base = illumina_R2_read_list[i].replace(".R2.fastq.gz", "")

            if (illumina_R1_base != illumina_R2_base):
                print("R1 and R2 entries do not match: \n  " + illumina_R1_read_list[i] + "\n  " + illumina_R2_read_list[i])
                sys.exit(1)
            read_base_list.append(os.path.split(illumina_R1_base)[1])
        return read_base_list

    def get_samples_with_illumina_PE_reads(self):
        return [sample for sample in self.get_sample_list() if len(self.get_illumina_PE_read_base_list(sample))]

    def remove_prefix_from_all_entries_in_list(self, _list, _prefix):
        result_list=[]
        for i in _list:
            if i.startswith(_prefix):
                # Remove the prefix using slicing
                result_list.append(i[len(_prefix):])
            else:
                # No match, keep the original string
                result_list.append(i)
        return result_list

    def get_illumina_PE_read_arguments(self, sample, argument_R1_prefix='', argument_R2_prefix=''):
        illumina_R1_read_lists = self.get_file_list(sample, "illumina-R1")
        illumina_R2_read_lists = self.get_file_list(sample, "illumina-R2")
        if len(illumina_R1_read_lists) != len(illumina_R2_read_lists):
            print("Sample does not have an equal number of R1 and R2 entries: " + sample)
            sys.exit(1)

        arg_list = []
        for i in range(0,len(illumina_R1_read_lists)):
            arg_list.append(argument_R1_prefix + illumina_R1_read_lists[i])
            arg_list.append(argument_R2_prefix + illumina_R2_read_lists[i])
        return " ".join(arg_list)

    #### The behavior of these functions depends on few factors
    #
    # MODE 1: argument prefix is 'references' - looks for files defined in data.csv (and uses their endings/format)
    # MODE 2: argument prefix is  <anything else>  - looks for files ending with the provided suffix
    #    If a suffix is not provided, then it will look specifically for files with endings in this order:
    #       *.gb, *.gbk, *.genbank, *.gff3, *.gff, *.fa, *.fna, *.fasta
    #    If none are found, it will default to 'fasta'
    # In the future, it would be nice to expand this fuzzy matching to allow different endings within a requested family

    preferred_reference_suffix_list = [
        'gb', 'gbk', 'genbank', 'gff3', 'gff', 'fa', 'fna'
    ]

    def add_preferred_reference_suffix(self, base_name):
        for prs in self.preferred_reference_suffix_list:
            test_name = base_name + "." + prs
            if (os.path.exists(test_name)):
                return(test_name)

        #Default if not found
        return base_name + ".fasta"

    def get_reference_arguments(self, sample, argument_prefix='', reference_suffix=None):
        if (self.reference_prefix == "references"):
            return " ".join([argument_prefix + os.path.join(self.reference_prefix, item) for item in self.get_file_list(sample, "reference")])
        else:
            if reference_suffix != None:
                return argument_prefix + os.path.join(self.reference_prefix, "{}.{}".format(sample, reference_suffix)) 
            else:
                return argument_prefix + self.add_preferred_reference_suffix(os.path.join(self.reference_prefix, sample))

    def get_specified_reference_list(self, reference_prefix, sample, reference_suffix=None):
        if (reference_prefix== "references"):
            return [ os.path.join(reference_prefix, item) for item in self.get_file_list(sample, "reference")]
        else:
            if reference_suffix != None:
                return os.path.join(reference_prefix, "{}.{}".format(sample, reference_suffix)) 
            else:
                return self.add_preferred_reference_suffix(os.path.join(reference_prefix, sample))

    def get_reference_list(self, sample, reference_suffix=None):
        return self.get_specified_reference_list(self.reference_prefix, sample, reference_suffix)

    def get_all_reference_list(self, reference_suffix=None):
        reference_arguments = dict()
        for s in sample_info.get_sample_list():
            for r in self.get_reference_list(s, reference_suffix):
                reference_arguments[r] = 1
        return reference_arguments.keys()

    def get_main_reference_file(self, reference_suffix=None):
        if reference_suffix != None:
            return "reference.{}".format(reference_suffix)
        else:
            return self.add_preferred_reference_suffix("reference")

    def get_local_path_list(self):
        return self.local_to_remote_path_mapping.keys()

    def get_remote_path_from_local_path(self, this_local_path):
        return self.local_to_remote_path_mapping[this_local_path]
        
    def get_sample_list(self):
        return list(self.sample_list)

    def set_reference_prefix(self, in_reference_prefix):
        self.reference_prefix = in_reference_prefix
    
    def get_reference_prefix(self):
        return self.reference_prefix

    def select_samples(self, sample_list):

        # Do this to avoid duplicates
        self.sample_set = set(sample_list)
        self.sample_list = list(self.sample_set)
        self.option_lists_by_sample_by_type = {key: self.option_lists_by_sample_by_type[key] for key in sample_list if key in self.option_lists_by_sample_by_type}
        self.file_lists_by_sample_by_type = {key: self.file_lists_by_sample_by_type[key] for key in sample_list if key in self.file_lists_by_sample_by_type}

        # Must also prune local_paths to avoid downloads
        local_path_set = set({})
        for this_sample in sample_list:
            if this_sample in self.local_paths_by_sample.keys():
                local_path_set.update(self.local_paths_by_sample[this_sample])

        self.local_to_remote_path_mapping = {key: self.local_to_remote_path_mapping[key] for key in local_path_set if key in self.local_to_remote_path_mapping}

#### Initialize from config values
brefito_config = {key.upper(): value for key, value in config.items()}

data_csv_name = 'data.csv'
if "data_csv" in config:
    data_csv_name = config["data_csv"]

print("\nSAMPLE INFORMATION\n")

sample_info = None

if sample_info == None:
    try:
        print("Attempting to load Sample Info from CSV file: " + data_csv_name)
        sample_info = SampleInfo(data_csv_name)
    except FileNotFoundError:
        print("  Could not find/load Sample Info file.")
    except Exception as e:
        print(f"Caught a general exception: {e}")
        import traceback
        traceback.print_exc()

if sample_info == None:
    print("Attempting to load Sample Info from directory structure: " + data_csv_name)
    sample_info = SampleInfo()

# Typically "references" or "assemblies" but can be overridden
if 'REFERENCES' in brefito_config.keys():
    sample_info.set_reference_prefix(brefito_config['REFERENCES'])

# Now set exactly which samples to use if option is provided
if 'SAMPLES' in brefito_config.keys():
    sample_info.select_samples(brefito_config['SAMPLES'].split("_,_"))

sample_info.print_file_lists()
print()