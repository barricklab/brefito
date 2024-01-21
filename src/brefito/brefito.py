#!/usr/bin/env python3

def main():

    import glob
    import os.path
    import argparse
    import re
    import subprocess


    # Where are the smk rules files?
    import brefito
    import importlib.resources
    brefito_package_path = importlib.resources.files(brefito)
    rules_path = os.path.join(brefito_package_path, "workflow", "rules")

    # What command did we choose
    parser = argparse.ArgumentParser(
                        prog='brefito',
                        description='wrapper for bacterial reference genome assembly, polishing, and annotation',
                        epilog='')

    parser.add_argument('-p', '--path', default='.', type=str)      # option that takes a value
    parser.add_argument('--config', action='append', default=[])         
    parser.add_argument('--rerun-incomplete', action='store_true') 
    parser.add_argument('--unlock', action='store_true') 
    parser.add_argument('--no-temp', action='store_true') 
    parser.add_argument('--keep-going', action='store_true') 
    parser.add_argument("--dry-run", '-r', action='store_true') 
    parser.add_argument("--resources", action='append', default=[]) 
    parser.add_argument('command', type=str)           # positional argument
    parser.add_argument('samples', nargs='?', type=str)           # positional argument

    args = parser.parse_args()
    base_path = args.path
    command_to_run = args.command.lower()
    samples_to_run = args.samples
    config_options_list = args.config
    resource_options_list = args.resources

    # What command did we run?
    print("Command: " + command_to_run)
    if (samples_to_run != None):
        print("Samples:" + samples_to_run)
    else:
        print("Samples: all")
    print("Base path: " + base_path)
    print("Config options: " + " ".join(config_options_list))
    print("Resource options: " + " ".join(resource_options_list))

### BEGIN MOVED

    # What files are available?

    def find_matching_input_files(in_base_path, in_file_ending):
        existing_files=glob.glob(os.path.join(in_base_path, "*."+in_file_ending))
        #print(os.path.join(base_path,"input", "*."+file_ending))
        matching_input_files = {}
        for this_input_file in existing_files:
            this_file_name=os.path.basename(this_input_file)
            #print(this_file_name)
            short_name = re.findall("(.+)\." + re.escape(in_file_ending), this_file_name)
            matching_input_files[short_name[0]] = this_input_file

        return(matching_input_files)

    assemblies_path = "assemblies"
    reference_assemblies_path = "references"
    sample_assemblies_path = "samples"
    nanopore_input_path = "nanopore_reads"
    illumina_input_path = "illumina_reads"

    print()
    print("Nanopore read files found (*.fastq.gz) in " + nanopore_input_path)
    print()

    input_nanopore_files = find_matching_input_files(nanopore_input_path, "fastq.gz")
    for k in input_nanopore_files: print("    " + k)
    if (len(input_nanopore_files.items()) == 0) : print("    " + "NONE FOUND")

    print()
    print("Paired-end Illumina read files found (*.R[1/2].gz) in " + illumina_input_path)
    print()

    input_illumina_1_files = find_matching_input_files(illumina_input_path, "R1.fastq.gz")
    input_illumina_2_files = find_matching_input_files(illumina_input_path, "R2.fastq.gz")

    for key in set( list(input_illumina_1_files.keys()) + list(input_illumina_2_files.keys()) ):
        if key in input_illumina_1_files.keys(): print ("    " + str(key) + " : " + input_illumina_1_files[key])
        if key in input_illumina_2_files.keys(): print ("    " + str(key) + " : " + input_illumina_2_files[key])

        # if (command_to_run != "download-reads-lftp") and (command_to_run != "download-data-lftp"):
        #     assert key in input_illumina_2_files.keys(), "Error: Matching R2 file does not exist"
        #     assert key in input_illumina_1_files.keys(), "Error: Matching R1 file does not exist"                                   

    if (len(input_illumina_1_files.items()) == 0) and (len(input_illumina_2_files.items()) == 0): 
        print("    " + "NONE FOUND")

    print()
    print("Genome assembly files found (*.fasta) in " + assemblies_path)
    print()
    input_assembly_files = find_matching_input_files(assemblies_path, "fasta")
    for (k, v) in input_assembly_files.items(): print("    " + k + " : " + v)
    if (len(input_assembly_files.items()) == 0) : print("    " + "NONE FOUND")

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
    print("Reference genome assembly files found (*.fasta) in " + reference_assemblies_path)
    print()
    input_reference_assembly_files = find_matching_input_files(reference_assemblies_path, "fasta")
    for (k, v) in input_reference_assembly_files.items(): print("    " + k + " : " + v)
    if (len(input_reference_assembly_files.items()) == 0) : print("    " + "NONE FOUND")


    print()
    print("Sample genome assembly files found (*.fasta) in " + sample_assemblies_path)
    print()
    input_sample_assembly_files = find_matching_input_files(sample_assemblies_path, "fasta")
    for (k, v) in input_sample_assembly_files.items(): print("    " + k + " : " + v)
    if (len(input_sample_assembly_files.items()) == 0) : print("    " + "NONE FOUND")

#### END MOVED

    snakemake_plus_common_options = ["snakemake", "--use-conda", "--cores", "all"]
    if args.rerun_incomplete:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--rerun-incomplete"]
    if args.unlock:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--unlock"]
    if args.no_temp:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--no-temp"]
    if args.keep_going:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--keep-going"]
    if args.dry_run:
        snakemake_plus_common_options = snakemake_plus_common_options + ["-r"]


    # What are appropriate targets for the command we are running?
    valid_command_found = False
    smk_targets = []

    if command_to_run == "polish-polypolish" or command_to_run == "polish-polca" or command_to_run == "polish-medaka":
        valid_command_found = True 
        smk_targets = [ d + ".polished" for d in input_assembly_files.values() ]

    # We use one smk with different targets for these three
    if command_to_run == "evaluate-breseq" or command_to_run == "evaluate-breseq-nanopore-reads":
        smk_targets = smk_targets + [ "evaluate/breseq_output/{}_nanopore_reads".format(key) for key in input_assembly_files.keys() ]
    if command_to_run == "evaluate-breseq" or command_to_run == "evaluate-breseq-illumina-reads":
        smk_targets = smk_targets + [ "evaluate/breseq_output/{}_illumina_reads".format(key) for key in input_assembly_files.keys() ]
    if command_to_run == "evaluate-breseq" or command_to_run == "evaluate-breseq-illumina-reads" or command_to_run == "evaluate-breseq-nanopore-reads":
        command_to_run = "evaluate-breseq"
        valid_command_found = True 

    # Set this globally
    resource_options_list = ["connections=1"] + resource_options_list

    if command_to_run == "download-reads-lftp":
        valid_command_found = True 

    if command_to_run == "download-data":
        valid_command_found = True 

    if command_to_run == "assemble-flye":
        smk_targets = []
        valid_command_found = True 

    if command_to_run == "assemble-unicycler":
        valid_command_found = True 

    if command_to_run == "assemble-unicycler-csv":
        valid_command_found = True 

    if command_to_run == "annotate-references":
        smk_targets = smk_targets + [ "output/annotated_references/{}.gbk".format(key) for key in input_reference_assembly_files.keys() ]
        valid_command_found = True 

    if command_to_run == "predict-mutations-breseq":

        #input_illumina_1_files.keys()
        #input_read_file_set = set(input_illumina_1_files.keys()) | set(input_illumina_2_files.keys()) | set(input_nanopore_files.keys())
        #smk_targets = smk_targets + [ "output/{}".format(key) for key in input_read_file_set ]
        valid_command_found = True 

    if command_to_run == "predict-mutations-breseq-csv":
       valid_command_found = True 

    if command_to_run == "evaluate-aligned-reads" or command_to_run == "align-for-igv":

        command_to_run = "evaluate-aligned-reads"
        # # Indexed FASTA
        # smk_targets = smk_targets + [ "evaluate/aligned_reads/{}/reference.fasta.fai".format(key) for key in input_assembly_files ]

        # #  Indexed nanopore BAM
        # for key in input_assembly_files:
        #     if (key in input_nanopore_files):
        #         smk_targets = smk_targets + [ "evaluate/aligned_reads/{}/nanopore.bam.bai".format(key) ]

        # #  Indexed illumina BAM
        # for key in input_assembly_files:
        #     if (key in input_illumina_1_files):
        #         smk_targets = smk_targets + [ "evaluate/aligned_reads/{}/illumina.bam.bai".format(key) ]

        valid_command_found = True 

    if command_to_run == "merge-references":
        valid_command_found = True 

    if command_to_run == "merge-trimmed-reads":
        valid_command_found = True 

    if command_to_run == "merge-reads":
        valid_command_found = True 


    # When we normalize, we need to know it will work = the same number of contigs as reference
    normalize_assembly_files = {}
    if command_to_run == "normalize-assemblies":
        from Bio import SeqIO
        assert input_main_reference_assembly_file != None, "Main reference assembly required for this command!" 
        input_main_reference_seqs = []
        for record in SeqIO.parse(input_main_reference_assembly_file, "fasta"):
            input_main_reference_seqs.append({'id' : record.id, 'seq' : record.seq})
        num_reference_contigs = len(input_main_reference_seqs)

        for input_assembly_file_key in input_assembly_files.keys():
            input_assembly_seqs = []
            input_assembly_file = input_assembly_files[input_assembly_file_key]
            for record in SeqIO.parse(input_assembly_file, "fasta"):
                input_assembly_seqs.append({'id' : record.id, 'seq' : record.seq})
            if len(input_assembly_seqs) == num_reference_contigs:
                smk_targets.append("assemblies/{}.fasta.normalized".format(input_assembly_file_key))
                normalize_assembly_files[input_assembly_file_key] = input_assembly_file
        valid_command_found = True 

    if command_to_run == "compare-syri":
        for s in input_sample_assembly_files.keys():
            for r in input_reference_assembly_files.keys():
                smk_targets = smk_targets + [ "comparisons/" + r + "/" + s + ".pdf" ]
        valid_command_found = True 

    if command_to_run == "compare-mummer":
        command_to_run = "compare-syri"
        for s in input_sample_assembly_files.keys():
            for r in input_reference_assembly_files.keys():
                smk_targets = smk_targets + [ "02_mummer_results/" + r + "/" + s + ".coords" ]
        valid_command_found = True 

    if command_to_run == "evaluate-nanopore-reads":
        smk_targets = smk_targets + [ "nanopore_read_stats/{}".format(key) for key in input_nanopore_files ]
        valid_command_found = True 

    if command_to_run == "evaluate-inspector":
        smk_targets = smk_targets + [ "inspector_assembly_evaluation/{}".format(key) for key in input_assembly_files ]
        valid_command_found = True 

    if command_to_run == "evaluate-coverage":
        smk_targets = smk_targets + [ "evaluate/coverage_plots/nanopore/{}".format(key) for key in input_assembly_files ]
        valid_command_found = True 

    if command_to_run == "evaluate-breseq-coverage":
        smk_targets = smk_targets + [ "evaluate/coverage_plots/breseq_nanopore/{}".format(key) for key in input_assembly_files ]
        valid_command_found = True 

    if command_to_run == "evaluate-isescan":
        smk_targets = smk_targets + [ "evaluate/isescan/{}.csv".format(key) for key in input_assembly_files ]
        valid_command_found = True 
    
    if command_to_run == "evaluate-soft-clipping":

        smk_targets = smk_targets + [ "evaluate/soft_clipping_summary/nanopore/{}_soft_clipping_summary.csv".format(key) for key in input_assembly_files ]
        config_options_list.append("brefito_package_path=" + str(brefito_package_path))
        valid_command_found = True 

    if command_to_run == "evaluate-redotable":
        smk_targets = smk_targets + [ "evaluate/dot_plot/{}.svg".format(key) for key in input_assembly_files ]
        #config_options_list.append("brefito_package_path=" + str(brefito_package_path))
        valid_command_found = True 


    #################################################
    ### trycycler trifecta
    #################################################

    if command_to_run == "trycycler-assemble":
        valid_command_found = True 
        smk_targets = [ "05_trycycler/" + d + "/done" for d in input_nanopore_files.keys() ]
        #resource_options_list = resource_options_list + ["necats=4"]

    if command_to_run == "trycycler-reconcile":
        valid_command_found = True 
        input_files=glob.glob("05_trycycler/*/cluster_*")
        for this_input_file in input_files:
            smk_targets.append(os.path.join(this_input_file, "2_all_seqs.fasta"))

    if command_to_run == "trycycler-consensus":
        valid_command_found = True 
        smk_targets = [ "assemblies/" + d + ".fasta" for d in input_nanopore_files.keys() ]

    #################################################
    assert valid_command_found, "Command not recognized: " + command_to_run + "\nRun brefito -h to see valid options!"

    smk_file_path = os.path.join(rules_path, command_to_run + ".smk")
    target_options = ["-s", smk_file_path]

    config_options = []
    if len(config_options_list) > 0:
        config_options =  ["--config"]
        for i in config_options_list:
            si = i.split('=')
            config_options = config_options + [si[0] + '="' + si[1] + '"']

    # This lets us replace defaults with user specified resource settings
    resource_options = []
    if len(resource_options_list) > 0:
        resource_options =  ["--resources"]
        resource_options_dict = {}
        for i in resource_options_list:
            si = i.split('=')
            resource_options_dict[si[0]]=si[1] 
        for k in resource_options_dict:
            resource_options = resource_options + [k + '=' + resource_options_dict[k]]

    command = snakemake_plus_common_options + target_options + smk_targets + config_options + resource_options

    print()
    print("RUNNING COMMAND")
    print(" ".join(command))

    subprocess.run(command)

    ## Cleanup

    def copy_and_rename_assemblies(in_input_assembly_files, in_ending_to_remove, in_ending_to_add):
        for a in in_input_assembly_files:
            
            # Check for new file
            if os.path.isfile(a + "." + in_ending_to_remove):

                # Check if we are the original file
                if not os.path.isfile(a + ".1.original"):
                    subprocess.run(["cp", a, a + ".1.original"])

                # Rename the new one and replace the main one so we can iterate
                i=1
                while len(glob.glob(a + "." + str(i) + ".*")) == 1:
                    i = i + 1

                subprocess.run(["cp", a + "." + in_ending_to_remove, a + "." + str(i) + "." + in_ending_to_add])
                subprocess.run(["mv", a + "." + in_ending_to_remove, a])


    if command_to_run == "polish-polypolish":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "polypolish")
    elif command_to_run == "polish-polca":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "polca")
    elif command_to_run == "polish-medaka":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "medaka")
    elif command_to_run == "normalize-assemblies":
        copy_and_rename_assemblies(normalize_assembly_files.values(), "normalized", "normalized")


if __name__ == "__main__":
    main()