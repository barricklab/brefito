# CONFIG OPTIONS
#
# bookmark = bookmark (REQUIRED)
# data_csv = path to CSV file (default = 'data.csv')
# remote_path = (default = '.')
#
# DATA_CSV
#
# This file must have the following columns:
#  1 sample = main sample name, must match for reads that will be used together.
#             May be a shell-style glob (containing *, ? or [...]) that applies the
#             row to every sample named exactly anywhere in the file (before or after).
#  2 type = either a file type: 'nanopore', 'nanopore-fast', 'nanopore-hac', 'nanopore-sup',
#           'pacbio', 'pacbio-clr', 'pacbio-hifi', 'illumina-SE', 'illumina-R1', 'illumina-R2',
#           'illumina-PE', 'reference' (see LONG_READ_TYPES below; 'nanopore-raw' and
#           'nanopore-hq' are accepted as synonyms, see LONG_READ_TYPE_SYNONYMS)
#           or the name of any --config option, like BRESEQ_OPTIONS (case-insensitive)
#  3 setting = for file types, a way of downloading the file:
#                 lftp@bookmark://path/to/file (if bookmark is present and refers to a location)
#                 wget://example.com/path/to/file
#              for options, the value of the setting (blank = true for on/off flags)
#
# PER-SAMPLE OPTIONS
#
#   *,BRESEQ_OPTIONS,--brief-html-output
#   REL1166A,BRESEQ_OPTIONS,--polymorphism-mode
#   REL*,BRESEQ_OPTIONS_ADD,--junction-only
#   REL1166A,NO_DEFAULT_BRESEQ_OPTIONS,
#
# For KEY: a row naming the sample exactly beats a glob row; among rows at the same
# level the later one wins. Every matching KEY_ADD row is appended after that, in file
# order. On the command line, --config KEY=... overrides data.csv (with a warning) and
# --config KEY_ADD=... (or brefito --config-add KEY=...) is put in front of it.
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
import fnmatch

# Run accessions of the three INSDC archives: SRR (NCBI SRA), ERR (ENA) and DRR (DDBJ). The SRA
# serves all three, so an sra:// setting may name any of them.
SRA_RUN_ACCESSION_PREFIXES = ("SRR", "ERR", "DRR")

def is_sra_run_accession(name):
    return name.startswith(SRA_RUN_ACCESSION_PREFIXES) and name[3:].isdigit()

# Sample names containing any of these characters in data.csv are globs
SAMPLE_PATTERN_CHARACTERS = set("*?[")

def is_sample_pattern(name):
    return any(c in SAMPLE_PATTERN_CHARACTERS for c in name)

# Long-read types. Each has its own read directories, <type>-reads/, <type>-reads-trimmed/,
# <type>-reads-filtered/, ..., so the nanopore paths are exactly what they always were. A sample
# may have reads of at most one long-read type; get_long_read_type() names it.
#
# The nanopore types are named for the Dorado basecaller model that called the reads (fast, hac,
# sup). The hac and sup settings assume R10.4.1 reads.
#
#  nanopore       nanopore reads of unstated accuracy: run with the nanopore-fast settings
#  nanopore-fast  noisier nanopore reads: Dorado/Guppy fast model, or older R9.4.1 reads
#                 (synonym: nanopore-raw)
#  nanopore-hac   Dorado high-accuracy (hac) model reads, R10.4.1 (synonym: nanopore-hq)
#  nanopore-sup   Dorado super-accuracy (sup) model reads, R10.4.1 (Q20+, under ~5% error)
#  pacbio         PacBio reads of unstated kind: run with the pacbio-clr settings
#  pacbio-clr     PacBio RS II/Sequel continuous long reads
#  pacbio-hifi    PacBio HiFi/CCS reads
#
# The unqualified types use the settings for the noisier reads, since those still work (if not
# optimally) on accurate reads, while settings for accurate reads give poor results on noisy ones.

# How each tool is told which kind of long reads it has.
#  platform      "ont" or "pacbio": ONT-only steps (porechop, medaka, breseq -x) need "ont"
#  minimap2      -x preset for mapping reads to a reference
#  minimap2_ava  -x preset for all-vs-all overlaps (miniasm)
#  flye          read-type option
#  canu          read-type option
#  autocycler    --read_type for autocycler helper
#  necat         NECAT can assemble these reads (it is built for nanopore reads)
#  metamdbg      metaMDBG can assemble these reads (it takes only HiFi and Q20+ nanopore reads)
#  medaka_model  default medaka model (MEDAKA_MODEL overrides it); None = medaka refuses them
_NANOPORE_FAST_SETTINGS = {"platform": "ont", "minimap2": "map-ont", "minimap2_ava": "ava-ont",
    "flye": "--nano-raw", "canu": "-nanopore", "autocycler": "ont_r9", "necat": True,
    "metamdbg": False, "medaka_model": "r941_min_sup_g507"}
# hac and sup reads get the same tool settings; only the medaka model (named for the
# basecaller model) differs
_NANOPORE_ACCURATE_SETTINGS = {"platform": "ont", "minimap2": "lr:hq", "minimap2_ava": "ava-ont",
    "flye": "--nano-hq", "canu": "-nanopore", "autocycler": "ont_r10", "necat": True,
    "metamdbg": True}
_PACBIO_CLR_SETTINGS = {"platform": "pacbio", "minimap2": "map-pb", "minimap2_ava": "ava-pb",
    "flye": "--pacbio-raw", "canu": "-pacbio", "autocycler": "pacbio_clr", "necat": False,
    "metamdbg": False, "medaka_model": None}

LONG_READ_SETTINGS = {
    "nanopore":      dict(_NANOPORE_FAST_SETTINGS, label="nanopore"),
    "nanopore-fast": dict(_NANOPORE_FAST_SETTINGS, label="nanopore (fast)"),
    "nanopore-hac":  dict(_NANOPORE_ACCURATE_SETTINGS, medaka_model="r1041_e82_400bps_hac_v4.3.0",
                          label="nanopore (high accuracy, hac)"),
    "nanopore-sup":  dict(_NANOPORE_ACCURATE_SETTINGS, medaka_model="r1041_e82_400bps_sup_v4.3.0",
                          label="nanopore (super accuracy, sup)"),
    "pacbio":       dict(_PACBIO_CLR_SETTINGS, label="PacBio"),
    "pacbio-clr":   dict(_PACBIO_CLR_SETTINGS, label="PacBio CLR"),
    "pacbio-hifi":  {"platform": "pacbio", "minimap2": "map-hifi", "minimap2_ava": "ava-pb",
                     "flye": "--pacbio-hifi", "canu": "-pacbio-hifi", "autocycler": "pacbio_hifi",
                     "necat": False, "metamdbg": True, "medaka_model": None,
                     "label": "PacBio HiFi"},
}

LONG_READ_TYPES = list(LONG_READ_SETTINGS)

# Older names accepted in data.csv. They are turned into the canonical name as data.csv is read,
# so only the names in LONG_READ_TYPES appear in directories, file names and wildcards.
LONG_READ_TYPE_SYNONYMS = {
    "nanopore-raw": "nanopore-fast",
    "nanopore-hq":  "nanopore-hac",
}

def canonical_long_read_type(name):
    """The LONG_READ_TYPES name for a long-read type or its synonym; None if it is neither."""
    name = LONG_READ_TYPE_SYNONYMS.get(name, name)
    return name if name in LONG_READ_SETTINGS else None

NANOPORE_TYPES = [t for t in LONG_READ_TYPES if LONG_READ_SETTINGS[t]["platform"] == "ont"]
PACBIO_TYPES = [t for t in LONG_READ_TYPES if LONG_READ_SETTINGS[t]["platform"] == "pacbio"]

# Wildcard regexes, longest names first so a prefix (nanopore) never cuts off a longer name
def _type_regex(names):
    return "|".join(sorted(names, key=len, reverse=True))

LONG_READ_TYPE_REGEX = _type_regex(LONG_READ_TYPES)
NANOPORE_TYPE_REGEX = _type_regex(NANOPORE_TYPES)
PACBIO_TYPE_REGEX = _type_regex(PACBIO_TYPES)

def long_read_setting(long_read_type, key):
    return LONG_READ_SETTINGS[long_read_type][key]

def is_nanopore_type(long_read_type):
    return long_read_type in NANOPORE_TYPES

# Aligned long reads are named <prefix>.<reads>.bam: nanopore_reads, nanopore_hac_reads, pacbio_clr_reads, ...
def long_read_bam_prefix(long_read_type):
    return long_read_type.replace("-", "_") + "_reads"

LONG_READ_BAM_PREFIXES = {long_read_bam_prefix(t): t for t in LONG_READ_TYPES}
LONG_READ_BAM_PREFIX_REGEX = _type_regex(LONG_READ_BAM_PREFIXES)

# --config keys brefito sets itself; data.csv may not set these
RESERVED_CONFIG_KEYS = set({"SAMPLES", "REFERENCES", "DATA_CSV", "BREFITO_PACKAGE_PATH", "QUERY_FILE_PATH", "BRESEQ_PRERELEASE", "BRESEQ_PATH"})

def rreplace(s, old, new, maxreplace):
    return(new.join(s.rsplit(old, maxreplace)))

def rreplace1(s, old, new):
    return(rreplace(s,old,new,1))

class SampleInfo():
    file_info = []

    sample_set = set({})
    sample_list = []

    reads_with_bad_endings = False

    # Dictionary of lists with keys that are types of data and then lists of files
    file_lists_by_sample_by_type = {}

    # data.csv rows that set options, in file order. Each is a dict with keys
    # pattern, is_pattern, key (upper case), value, row_index. Matched lazily
    # against whatever sample is being looked up.
    option_rows = []

    # File rows whose sample is a glob. Expanded once the whole CSV is read.
    pending_pattern_file_rows = []

    # Globs (file or option rows) that matched no sample
    unmatched_sample_patterns = []

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
            first_column_name = data_reader.fieldnames[0] if data_reader.fieldnames else None
            for row_index, row in enumerate(data_reader):
                #print(row)

                first_column_value = row.get(first_column_name)

                # Skip blank lines
                if first_column_value is None or first_column_value.strip() == "":
                    continue

                #Ignore commented rows, where first item begins with #
                if first_column_value.lstrip().startswith("#"):
                    continue

                # Add 'path' for 'setting' for backwards compatibility
                if 'path' in row:
                    if ('setting' in row):
                        print("Both 'path' and 'setting' found for row. Only using 'setting'.")
                    else:
                        row["setting"] = row["path"]

                # A row that did not split into the expected columns leaves 'type'
                # and/or 'setting' as None. Fail with a clear, actionable message
                # instead of a cryptic error deeper in.
                if row.get('type') is None or row.get('setting') is None:
                    hint = ""
                    if "\t" in first_column_value:
                        hint = (" This row looks TAB-separated, but brefito data files must be "
                                "COMMA-separated (the header and every data row). Re-save it "
                                "with commas between the sample, type, and setting columns.")
                    raise ValueError(
                        "could not parse a row into the required 'sample,type,setting' "
                        "columns: " + repr(first_column_value) + "." + hint)

                # Copy over certain paths to allow for some choices.
                file_name = os.path.basename(row['setting'])
                base_file_name = os.path.basename(row['setting'])

                # A synonym (nanopore-hq, ...) is stored under its canonical name
                long_read_type = canonical_long_read_type(row['type'])
                if long_read_type is not None:
                    simplified_read_file_base_name = self.get_simplified_read_file_base_name(base_file_name)
                    local_path = os.path.join(long_read_type + "-reads", simplified_read_file_base_name + ".fastq.gz")
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : long_read_type,
                        'file_type' : long_read_type,
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina") or (row['type'] == "illumina-SE"):
                    simplified_read_file_base_name = self.get_simplified_read_file_base_name(base_file_name)
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".SE.fastq.gz")
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "illumina-SE",
                        'file_type' : "illumina",
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R1"):
                    simplified_read_file_base_name = self.get_simplified_read_file_base_name(base_file_name)
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".R1.fastq.gz")
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "illumina-R1",
                        'file_type' : "illumina",
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                elif (row['type'] == "illumina-R2"):
                    simplified_read_file_base_name = self.get_simplified_read_file_base_name(base_file_name)
                    local_path = os.path.join("illumina-reads", simplified_read_file_base_name + ".R2.fastq.gz")
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "illumina-R2",
                        'file_type' : "illumina",
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })

                # paired end reads
                elif row['type'] in ["illumina-pair", "illumina-paired", "illumina-PE"]:
                    
                    # Don't do the automatic read name simplification step for sra:// run accessions
                    # (SRR, ERR or DRR)
                    if row['setting'].startswith("sra://") and is_sra_run_accession(row['setting'][6:]):

                        #Make sure we cut off the sra:// part using slice
                        simplified_read_file_base_name_1 = row['setting'][6:] # this value is already simplified
                        simplified_read_file_base_name_2 = row['setting'][6:] # this value is already simplified

                        remote_path_1 = row['setting']
                        remote_path_2 = row['setting']

                    # Normal case non-SRA
                    else:
                        if not "{1|2}" in row['setting']:
                            print("Did not find {1|2} in name of illumina-PE entry:" + row['setting'])
                            print("Skipping this entry.")
                            continue
                                                
                        simplified_read_file_base_name_1 = self.get_simplified_read_file_base_name(base_file_name.replace("{1|2}", "1"))
                        simplified_read_file_base_name_2 = self.get_simplified_read_file_base_name(base_file_name.replace("{1|2}", "2"))

                        remote_path_1 = row['setting'].replace("{1|2}", "1")
                        remote_path_2 = row['setting'].replace("{1|2}", "2")

                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "illumina-R1",
                        'file_type' : "illumina",
                        'remote_path' : remote_path_1,
                        'local_path' : os.path.join("illumina-reads", simplified_read_file_base_name_1 + ".R1.fastq.gz")
                        })
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "illumina-R2",
                        'file_type' : "illumina",
                        'remote_path' : remote_path_2,
                        'local_path' : os.path.join("illumina-reads", simplified_read_file_base_name_2 + ".R2.fastq.gz")
                        })

                elif (row['type'] == "reference"):
                    local_path = os.path.join("references", file_name)
                    self.add_csv_file({
                        'sample' : row['sample'],
                        'type' : "reference",
                        'file_type' : "reference",
                        'remote_path' : row['setting'],
                        'local_path' : local_path
                        })
                else:
                    # Any other type is the name of a --config option
                    option_key = row['type'].strip().upper()
                    option_base_key = option_key[:-len("_ADD")] if option_key.endswith("_ADD") else option_key
                    if option_base_key in RESERVED_CONFIG_KEYS:
                        raise ValueError("option '" + row['type'] + "' cannot be set in a data file "
                                         "(it is set by brefito itself). Use the brefito command line instead.")
                    self.option_rows.append({
                        'pattern' : row['sample'].strip(),
                        'is_pattern' : is_sample_pattern(row['sample'].strip()),
                        'key' : option_key,
                        'value' : row['setting'].strip(),
                        'row_index' : row_index,
                        })

                #else:
                #    print("Skipping row with unknown type:" + row['type'])

        # Globs apply to every sample named exactly anywhere in the file, so they
        # can only be expanded now that the whole file has been read.
        exact_samples = sorted(self.sample_set)
        for row in self.pending_pattern_file_rows:
            matched_samples = [s for s in exact_samples if fnmatch.fnmatchcase(s, row['sample'])]
            if len(matched_samples) == 0:
                self.unmatched_sample_patterns.append(row['sample'] + " (" + row['type'] + " " + str(row['remote_path']) + ")")
            for this_sample in matched_samples:
                this_row = copy.deepcopy(row)
                this_row['sample'] = this_sample
                self.add_file(this_row)

        for option_row in self.option_rows:
            if option_row['is_pattern'] and not any(fnmatch.fnmatchcase(s, option_row['pattern']) for s in exact_samples):
                self.unmatched_sample_patterns.append(option_row['pattern'] + " (" + option_row['key'] + ")")

        self.sample_list = list(self.sample_set)
        return True

    # File rows from data.csv: glob samples are held back until the whole file is read
    def add_csv_file(self, row):
        if is_sample_pattern(row['sample']):
            self.pending_pattern_file_rows.append(row)
        else:
            self.add_file(row)

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
        illumina_input_path = "illumina-reads"

        # Only canonical type names are directory names. Reads in a directory named for a
        # synonym (nanopore-hq-reads/, ...) would be silently left out, so stop instead.
        for synonym, long_read_type in LONG_READ_TYPE_SYNONYMS.items():
            if len(self.find_matching_input_files(synonym + "-reads", ".fastq.gz")):
                print("ERROR: Found read files in " + synonym + "-reads/, but '" + synonym + "' reads are now")
                print("  typed '" + long_read_type + "' and are looked for in " + long_read_type + "-reads/.")
                print("  Rename the directory to " + long_read_type + "-reads/ (or list the files in a data.csv).")
                exit(1)

        for long_read_type in LONG_READ_TYPES:
            long_read_input_path = long_read_type + "-reads"
            input_long_read_files = self.find_matching_input_files(long_read_input_path, ".fastq.gz")
            # Only nanopore-reads/ is reported when empty, as it always was; the other
            # long-read directories are mentioned only when they hold reads.
            if long_read_type != "nanopore" and len(input_long_read_files) == 0:
                continue
            print()
            print(long_read_setting(long_read_type, "label") + " read files found (*.fastq.gz) in " + long_read_input_path)
            print()
            for k in input_long_read_files: print("    " + k)
            if (len(input_long_read_files.items()) == 0) : print("    " + "NONE FOUND")
            for k, i in input_long_read_files.items():
                self.add_file({
                            'sample' : k,
                            'type' : long_read_type,
                            'file_type' : long_read_type,
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

        # add_file() populates self.sample_set; mirror init_from_CSV_file and
        # publish it as sample_list so get_sample_list() works under autodetect.
        self.sample_list = list(self.sample_set)

    def add_file(self, row):

        # Do checks that lead to not adding to the list first
        #if row['remote_path']!= None and self.remote_path_is_duplicate(row['remote_path']):
        #    print("Duplicate remote path encountered. Skipping entry:")
        #    print(row)
        #    return

        # Glob samples from data.csv are expanded before reaching here (see add_csv_file)
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

    # Entries that name the same source share one local file. It is downloaded
    # once and every sample that lists it uses that same copy. Only genuinely
    # different sources that would land on the same local name get a "-1", "-2",
    # ... suffix added to keep them separate.
    def deconflict_paths(self, this_row):

        candidate_path = this_row['local_path']
        file_name, file_extension = self.split_filename_and_extension(candidate_path)

        i = 0
        while True:
            claimed_path = self.find_claimed_local_path(candidate_path)

            # Nobody is using this name yet
            if claimed_path == None:
                this_row['local_path'] = candidate_path
                return

            # Same source file - reuse the local copy already claimed for it
            if self.local_to_remote_path_mapping[claimed_path] == this_row['remote_path']:
                this_row['local_path'] = claimed_path
                return

            # Different sources - rename this one so they don't overwrite each other
            i = i + 1
            new_candidate_path = file_name + "-" + str(i) + file_extension
            print("WARNING: Two different sources map to the same local file name:")
            print("    " + str(self.local_to_remote_path_mapping[claimed_path]) + " => " + claimed_path)
            print("    " + str(this_row['remote_path']) + " => " + candidate_path)
            print("  Renaming the second one to: " + new_candidate_path)
            candidate_path = new_candidate_path

    # local_to_remote_path_mapping is keyed by the exact local path, but names that
    # differ only in case would collide on a case-insensitive file system.
    def find_claimed_local_path(self, local_path):
        for claimed_path in self.local_to_remote_path_mapping.keys():
            if claimed_path.lower() == local_path.lower():
                return claimed_path
        return None

    # A bare SRA run accession (what an sra:// entry gives us) has no .fastq.gz ending to
    # strip and is already in simplified form. A real file whose name merely begins with an
    # accession (SRR13512703_1.fastq.gz) must still be simplified, or the two files of a pair
    # end up with different base names and no paired rule can be built for them.
    def is_bare_sra_run_accession(self, in_read_name):
        return is_sra_run_accession(in_read_name)

    ## We want the read names to be standardized... this should do it in most cases
    def get_simplified_read_file_base_name(self, in_read_name):

        if self.is_bare_sra_run_accession(in_read_name):
            return in_read_name

        # Valid read file names must end in *.fastq.gz for workflows to function!
        no_ending_read_name = rreplace1(in_read_name, ".fastq.gz", "")
        if no_ending_read_name == in_read_name:
            self.reads_with_bad_endings = True
            return "INVALID---[" + in_read_name + "]---"

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
            # A sample may have no associated files (e.g. samples requested at the CLI
            # when no data.csv is present and none were found by the directory scan).
            file_lists_by_type = self.file_lists_by_sample_by_type.get(this_sample, {})
            for this_type in sorted(file_lists_by_type.keys()):
                for this_item in file_lists_by_type[this_type]:
                    print("  " + this_type + " : " + this_item)
        #print(self.local_to_remote_path_mapping)

    # Nanopore reads of any accuracy (nanopore, nanopore-fast, nanopore-hac, nanopore-sup)
    def get_nanopore_read_list(self, sample):
        return [f for t in NANOPORE_TYPES for f in self.get_file_list(sample, t)]

    def get_nanopore_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.get_nanopore_read_list(sample)])

    def get_nanopore_read_base_list(self, sample):
        return [os.path.split(n.replace(".fastq.gz", ""))[1] for n in self.get_nanopore_read_list(sample)]

    def get_samples_with_nanopore_reads(self):
        return [sample for sample in self.get_sample_list() if len(self.get_nanopore_read_list(sample))]

    ## Long reads of any type (LONG_READ_TYPES). A sample has at most one type.

    def get_long_read_type(self, sample):
        """The long-read type of this sample's reads, or None if it has none."""
        types = [t for t in LONG_READ_TYPES if len(self.get_file_list(sample, t))]
        if len(types) > 1:
            raise ValueError("sample " + sample + " has reads of more than one long-read type (" +
                             ", ".join(types) + "); give each kind its own sample")
        return types[0] if types else None

    def get_long_read_list(self, sample):
        """This sample's long-read files (base names), whatever their type."""
        t = self.get_long_read_type(sample)
        return self.get_file_list(sample, t) if t else []

    def get_long_read_base_list(self, sample):
        return [os.path.split(n.replace(".fastq.gz", ""))[1] for n in self.get_long_read_list(sample)]

    def get_long_read_arguments(self, sample, argument_prefix=''):
        return " ".join([argument_prefix + item for item in self.get_long_read_list(sample)])

    def get_samples_with_long_reads(self, long_read_type=None):
        """Samples with long reads (of the given type, if one is given)."""
        return [s for s in self.get_sample_list() if self.get_long_read_type(s) is not None
                and (long_read_type is None or self.get_long_read_type(s) == long_read_type)]

    def get_long_read_path(self, sample, stage, suffix=".fastq"):
        """Per-sample long-read file at a processing stage, e.g. ('filtered', '.fastq') gives
        nanopore-reads-filtered/<sample>.fastq for a nanopore sample."""
        t = self.get_long_read_type(sample)
        if t is None:
            return None
        return t + "-reads-" + stage + "/" + sample + suffix

    def check_long_read_types(self):
        """Fail early, with a clear message, if a sample mixes long-read types."""
        for s in self.get_sample_list():
            self.get_long_read_type(s)

    def get_illumina_read_list(self, sample):
        return self.get_file_list(sample, "illumina-SE") + self.get_file_list(sample, "illumina-R1") + self.get_file_list(sample, "illumina-R2")

    def get_illumina_SE_read_list(self, sample):
        return self.get_file_list(sample, "illumina-SE")

    def get_illumina_SE_read_base_list(self, sample):
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
                print("  Names must match after the .R1/.R2 part is removed. Check the 'setting' column in data.csv.")
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
        self.file_lists_by_sample_by_type = {key: self.file_lists_by_sample_by_type[key] for key in sample_list if key in self.file_lists_by_sample_by_type}

        # Must also prune local_paths to avoid downloads
        local_path_set = set({})
        for this_sample in sample_list:
            if this_sample in self.local_paths_by_sample.keys():
                local_path_set.update(self.local_paths_by_sample[this_sample])

        self.local_to_remote_path_mapping = {key: self.local_to_remote_path_mapping[key] for key in local_path_set if key in self.local_to_remote_path_mapping}


    ## Per-sample options from data.csv

    def get_matching_option_rows(self, sample, key):
        return [r for r in self.option_rows if r['key'] == key and
                (fnmatch.fnmatchcase(sample, r['pattern']) if r['is_pattern'] else r['pattern'] == sample)]

    # Returns (replace_value_or_None, [add_values]) for an option key (upper case,
    # without _ADD). For KEY rows, one naming the sample exactly beats a glob; among
    # rows at the same level, the later row wins. All KEY_ADD rows apply, in file order.
    def get_csv_option(self, sample, key):
        key = key.upper()
        replace_rows = self.get_matching_option_rows(sample, key)
        replace_value = None
        if len(replace_rows) > 0:
            exact_rows = [r for r in replace_rows if not r['is_pattern']]
            winning_rows = exact_rows if len(exact_rows) > 0 else replace_rows
            replace_value = winning_rows[-1]['value']
        add_values = [r['value'] for r in self.get_matching_option_rows(sample, key + "_ADD")]
        return (replace_value, add_values)

    def csv_sets_option(self, sample, key):
        replace_value, add_values = self.get_csv_option(sample, key)
        return replace_value != None or len(add_values) > 0

    # Base option keys (no _ADD suffix) that data.csv sets for a sample
    def get_csv_option_keys(self, sample):
        keys = []
        for r in self.option_rows:
            base_key = r['key'][:-len("_ADD")] if r['key'].endswith("_ADD") else r['key']
            if base_key not in keys and self.csv_sets_option(sample, base_key):
                keys.append(base_key)
        return keys

    # Which samples list this read file (by base name or local path)?
    def get_samples_for_read_file(self, read_file_name):
        read_file_name = os.path.basename(read_file_name)
        samples = []
        for this_sample, files_by_type in self.file_lists_by_sample_by_type.items():
            for this_type, file_list in files_by_type.items():
                if read_file_name in file_list and this_sample not in samples:
                    samples.append(this_sample)
        return samples

#### Initialize from config values
brefito_config = {key.upper(): value for key, value in config.items()}

def join_option_values(values):
    return " ".join([str(v) for v in values if v != None and str(v).strip() != ""])

def get_config_value(key, default=None, sample=None):
    """Resolve a --config option, optionally for one sample.

    Order: the workflow default; replaced by data.csv KEY for this sample (exact name
    beats glob, later row beats earlier); then data.csv KEY_ADD values appended. A
    command-line --config KEY replaces all of that. A command-line KEY_ADD (also set by
    brefito --config-add KEY=...) is put in front. Returns default unchanged if nothing
    sets the key.
    """
    key = key.upper()
    value = default
    if sample != None and sample_info != None:
        csv_replace_value, csv_add_values = sample_info.get_csv_option(sample, key)
        if csv_replace_value != None:
            value = csv_replace_value
        if len(csv_add_values) > 0:
            value = join_option_values([value] + csv_add_values)
    if key in brefito_config:
        value = brefito_config[key]
    if key + "_ADD" in brefito_config:
        value = join_option_values([brefito_config[key + "_ADD"], value])
    return value

def get_uniform_config_value(key, default=None, samples=None):
    """Resolve an option that cannot vary between the given samples (all selected
    samples by default), e.g. one that picks a rule order or names output files. Exits
    with an error if data.csv gives them different values."""
    if samples == None:
        samples = sample_info.get_sample_list() if sample_info != None else []
    if len(samples) == 0:
        return get_config_value(key, default)
    values_by_sample = {s : get_config_value(key, default, s) for s in sorted(samples)}
    distinct_values = []
    for v in values_by_sample.values():
        if v not in distinct_values:
            distinct_values.append(v)
    if len(distinct_values) > 1:
        print()
        print("ERROR: Option " + key.upper() + " must be the same for these samples, but data.csv sets it differently:")
        for s, v in values_by_sample.items():
            print("    " + s + " : " + str(v))
        print("  Set it the same way for all of them in data.csv (e.g. with a '*' row),")
        print("  or override it for all samples with --config " + key.upper() + "=<value>.")
        exit(1)
    return distinct_values[0]

def value_is_true(v):
    if isinstance(v, bool):
        return v
    # A blank data.csv setting means "on" (the flag was listed without a value)
    if str(v).strip() == "":
        return True
    return str(v).strip().lower() not in ("0", "false", "no", "off", "f", "none")

def config_is_true(name, default=False, sample=None):
    """Interpret an on/off option as a boolean, consistently across workflows.

    Handles `--config NAME` (the brefito wrapper turns a value-less flag into "1"),
    `--config NAME=1|true|yes|on`, and `NAME=0|false|no|off`, plus the same values (or a
    blank setting, meaning true) in data.csv for `sample`. An unset option returns
    `default`. Use ONLY for genuine flags - not for value-bearing options.
    """
    v = get_config_value(name, None, sample)
    if v == None:
        return default
    return value_is_true(v)

def config_is_true_for_all_samples(name, default=False):
    """Like config_is_true for a flag that cannot vary between samples."""
    v = get_uniform_config_value(name, None)
    if v == None:
        return default
    return value_is_true(v)

BRESEQ_ENV = "../envs/breseq-prerelease.yml" if brefito_config.get("BRESEQ_PRERELEASE") else "../envs/breseq.yml"

# --breseq-path DIR: prefer breseq/gdtools from a local build over the conda package.
# Applied to every rule, so a locally built binary is used wherever one is called.
if brefito_config.get("BRESEQ_PATH"):
    print("Using breseq/gdtools from: " + brefito_config.get("BRESEQ_PATH"))
    shell.prefix('export PATH="' + brefito_config.get("BRESEQ_PATH") + ':$PATH"; ')

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
        # No such file: legitimately fall back to autodetecting from the directory
        # structure below (e.g. the default 'data.csv' simply isn't present).
        print("  Could not find Sample Info file '" + data_csv_name + "'.")
    except Exception as e:
        # The file exists but could not be parsed. Do NOT silently fall back to a
        # directory scan (that hides the real cause and yields a misleading "none
        # found" message) - surface the error clearly and exit.
        print("ERROR: Could not parse sample information file '" + data_csv_name + "':")
        print("  " + str(e))
        exit(1)

if sample_info == None:
    print("Attempting to load Sample Info from directory structure")
    sample_info = SampleInfo()

# Typically "references" or "assemblies" but can be overridden
if 'REFERENCES' in brefito_config.keys():
    sample_info.set_reference_prefix(brefito_config['REFERENCES'])

# Now set exactly which samples to use if option is provided
if 'SAMPLES' in brefito_config.keys():
    sample_info.select_samples(brefito_config['SAMPLES'].split("_,_"))

sample_info.print_file_lists()
print()

try:
    sample_info.check_long_read_types()
except ValueError as e:
    print("ERROR: " + str(e))
    exit(1)

def print_per_sample_config():
    selected_samples = sorted(sample_info.get_sample_list())
    samples_with_options = [s for s in selected_samples if len(sample_info.get_csv_option_keys(s)) > 0]
    if len(samples_with_options) > 0:
        print("PER-SAMPLE CONFIG (data.csv)")
        for this_sample in samples_with_options:
            print()
            print("SAMPLE :: " + this_sample)
            for this_key in sample_info.get_csv_option_keys(this_sample):
                this_value = str(get_config_value(this_key, "", this_sample))
                print("  " + this_key + " : " + (this_value if this_value.strip() != "" else "(blank = on)"))
        print()

    # Command-line --config KEY replaces data.csv settings for KEY: say so once per key
    overridden = {}
    for this_sample in selected_samples:
        for this_key in sample_info.get_csv_option_keys(this_sample):
            if this_key in brefito_config:
                overridden.setdefault(this_key, []).append(this_sample)
    for this_key, these_samples in overridden.items():
        print("WARNING: --config " + this_key + "=" + str(brefito_config[this_key]) + " on the command line overrides")
        print("  the data.csv setting of " + this_key + " for samples: " + ", ".join(these_samples))
        print("  Use --config-add " + this_key + "=... (or --config " + this_key + "_ADD=...) to add to it instead.")
        print()

    for this_pattern in sample_info.unmatched_sample_patterns:
        print("WARNING: data.csv sample pattern matched no samples: " + this_pattern)
    if len(sample_info.unmatched_sample_patterns) > 0:
        print()

print_per_sample_config()

if (sample_info.reads_with_bad_endings):
  print("ERROR: Read file names not of the expected form (*.fastq.gz) found.")
  print("Check above for offending file names in the form \"INVALID---[in_read_name]---*.fastq.gz\"")
  exit(1)

# A workflow may set SAMPLE_INFO_REQUIRED = False before including this file if it
# derives its samples elsewhere (see curate-ltee-clones). Default: required.
try: SAMPLE_INFO_REQUIRED
except NameError: SAMPLE_INFO_REQUIRED = True

if SAMPLE_INFO_REQUIRED and not sample_info.get_sample_list():
    print("ERROR: This workflow requires sample information, but none was found.")
    print("  Provide a data.csv (or a file via -d) describing your samples, or place input")
    print("  files in the expected directories (" + ", ".join(t + "-reads/" for t in LONG_READ_TYPES) + ",")
    print("  illumina-reads/, assemblies/, references/).")
    exit(1)
