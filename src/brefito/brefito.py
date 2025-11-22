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

    # Find all of the *.smk files in rules_path
    # These are valid by default
    valid_workflows = glob.glob(os.path.join(rules_path, "*.smk"))
    valid_workflows = sorted([os.path.basename(f).replace(".smk", "") for f in valid_workflows])

    valid_workflow_help = "Valid workflows are:\n  " + "\n  ".join(valid_workflows)

    # OLD => OLD should be deleted when migration complete
    assemblies_path = "assemblies"
    reference_assemblies_path = "references"
    sample_assemblies_path = "samples"
    nanopore_input_path = "nanopore_reads"
    illumina_input_path = "illumina_reads"
    # <= OLD should be deleted when migration complete

    # What command did we choose
    parser = argparse.ArgumentParser(
                        prog='brefito',
                        description='Wrapper script for Snakemake bacterial reference genome assembly, polishing, and annotation workflows.',
                        epilog=valid_workflow_help,
                        formatter_class=argparse.RawDescriptionHelpFormatter
                        )

    parser.add_argument('-p', '--path', default='.', type=str)

    #An additional way to specify these
    parser.add_argument('-r', '--references', type=str, help='Path to use for reference files. Default: references')

    #Snakemake passthroughs
    parser.add_argument('--cores', type=int, default=0, help='--cores argument passed through to Snakemake (0 = all)') # 0 means "all"
    parser.add_argument('--config', action='append', default=[], help='--config argument passed through to Snakemake. Individual workflows support different settings.')
    parser.add_argument('--resources', action='append', default=[], help='--resources argument passed through to Snakemake. Individual workflows support different settings.')
    parser.add_argument('--rerun-incomplete', action='store_true', help='argument passed through to Snakemake') 
    parser.add_argument('--rerun-triggers', action='append', default=['code', 'input', 'params'], help='argument passed through to Snakemake.') 
    parser.add_argument('--unlock', action='store_true', help='argument passed through to Snakemake') 
    parser.add_argument('--notemp', action='store_true', help='argument passed through to Snakemake')
    parser.add_argument('--keep-going', action='store_true', help='argument passed through to Snakemake') 
    parser.add_argument('--dry-run', action='store_true', help='argument passed through to Snakemake')

    # REQUIRED positional argument
    parser.add_argument('workflow', type=str)

     # OPTIONAL positional arguments
    parser.add_argument('samples', nargs='*', type=str)

    args = parser.parse_args()

    base_path = args.path
    workflow_to_run = args.workflow.lower()
    samples_to_run = args.samples

    config_options_list = args.config
    resource_options_list = args.resources
    references_argument = args.references

    # Is the workflow valid?

    # First check for workflows that take -* wildcards specifying 
    # the references and update the workflow and references accordingly

    # List longer overlapping matches first, so they are preferred
    # For example,  "align-reads-merged", then "align-reads".
    match = check_command_list_with_references(
        workflow_to_run, [
            "predict-mutations-breseq",
            "predict-mutations-minimap2-breseq",
            "evaluate-coverage-breseq",
            "predict-cnv-breseq",
            "compare-mutations-breseq",
            "coverage-plots-breseq",
            "align-reads",
            "check-soft-clipping",
            "mutate-genomes-gdtools",
            "annotate-genomes",
            "search-blast",
            "tabulate-ssrs-breseq"
            ]
    )

    if match['matched']:

        if match['references'] != None:
            if references_argument != None:
                print()
                print("OPTIONS WARNING")
                print("  Workflow suffix specified reference: " + match['references'])
                print("  Overrides command line option references " + references_argument)
            references_argument = match['references']

        workflow_to_run = match['workflow_to_run']


    # Now check whether it is valid
    assert workflow_to_run in valid_workflows, "Workflow not recognized: " + workflow_to_run + "\n" + valid_workflow_help

    # For search-blast-* Take the last unnamed argument, and there has to be one!
    # Use it to create the input file if it is not a file name
    search_query_option = ""
    if (workflow_to_run=="search-blast"):
        if len(samples_to_run)==0:
            print("For the search-blast-* command, you must provide at least one unnamed argument, which is either a DNA sequence or a FASTA file of DNA sequences.")
            return
        
        search_query_option = samples_to_run.pop(0)

        output_path = "blast"
        if references_argument:
            output_path += "-" + references_argument
        else:
            output_path += "-references"

        # User provided a raw sequence. Create a file with the default name
        if not os.path.exists(search_query_option):
            os.makedirs(output_path, exist_ok=True)
            default_filename = os.path.join(output_path, "query_sequence.fasta")
            with open(default_filename, 'w') as sequence_file:
                sequence_file.write(">input\n")
                sequence_file.write(search_query_option + "\n")
                search_query_option = default_filename
        
        # Delete old output files
        for s in samples_to_run:
            subprocess.run(["rm", os.path.join(output_path, "blast_" + s + "_output.html")])


    # Print out some details to help users debug bad command lines
    print("Workflow: " + workflow_to_run)
    if (samples_to_run != None):
        print("Samples:" + str(samples_to_run))
    else:
        print("Samples: all")
    print("Base path: " + base_path)
    print("Config options:")
    for i in config_options_list:
        print("  " + i)
    print("Resource options:")
    for i in resource_options_list:
        print("  " + i)



### BEGIN TODELETE WHEN MIGRATION COMPLETE

    # What files are available?

    def find_matching_input_files(in_base_path, in_file_ending):
        existing_files=glob.glob(os.path.join(in_base_path, "*."+in_file_ending))
        #print(os.path.join(base_path,"input", "*."+file_ending))
        matching_input_files = {}
        for this_input_file in existing_files:
            this_file_name=os.path.basename(this_input_file)
            #print(this_file_name)
            short_name = re.findall(r'(.+)\.' + re.escape(in_file_ending), this_file_name)
            matching_input_files[short_name[0]] = this_input_file

        return(matching_input_files)

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

        # if (workflow_to_run != "download-reads-lftp") and (workflow_to_run != "download-data-lftp"):
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

#### END TODELETE WHEN MIGRATION COMPLETE

    snakemake_plus_common_options = ["snakemake", "--use-conda"]
    if args.cores == 0:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--cores", "all"]
    else:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--cores", str(args.cores)]

    for rt in args.rerun_triggers:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--rerun-triggers", rt]

    if args.rerun_incomplete:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--rerun-incomplete"]
    if args.unlock:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--unlock"]
    if args.notemp:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--notemp"]
    if args.keep_going:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--keep-going"]
    if args.dry_run:
        snakemake_plus_common_options = snakemake_plus_common_options + ["-r"]

    # What are appropriate targets for the workflow we are running?
    smk_targets = []

    if workflow_to_run in ["polish-breseq", "polish-polypolish", "polish-polca", "polish-medaka"]:
        smk_targets = [ d + ".polished" for d in input_assembly_files.values() ]

    # Set this globally, putting it first means it can be overridden
    resource_options_list = ["connections=1"] + resource_options_list

    if workflow_to_run == "check-soft-clipping":
        config_options_list.append("brefito_package_path=" + str(brefito_package_path))

    # If not specified at command line or in workflow, set to default        
    if references_argument == None:
        references_argument = 'references';


    ## Commands that haven't yet been updated below --->

    # When we normalize, we need to know it will work = the same number of contigs as reference
    # normalize_assembly_files = {}
    # if workflow_to_run == "normalize-genomes":
    #     from Bio import SeqIO
    #     assert input_main_reference_assembly_file != None, "Main reference assembly required for this command!" 
    #     input_main_reference_seqs = []
    #     for record in SeqIO.parse(input_main_reference_assembly_file, "fasta"):
    #         input_main_reference_seqs.append({'id' : record.id, 'seq' : record.seq})
    #     num_reference_contigs = len(input_main_reference_seqs)

    #     for input_assembly_file_key in input_assembly_files.keys():
    #         input_assembly_seqs = []
    #         input_assembly_file = input_assembly_files[input_assembly_file_key]
    #         for record in SeqIO.parse(input_assembly_file, "fasta"):
    #             input_assembly_seqs.append({'id' : record.id, 'seq' : record.seq})
    #         if len(input_assembly_seqs) == num_reference_contigs:
    #             smk_targets.append("assemblies/{}.fasta.normalized".format(input_assembly_file_key))
    #             normalize_assembly_files[input_assembly_file_key] = input_assembly_file
    #     valid_command_found = True 

    if workflow_to_run == "compare-syri":
        for s in input_sample_assembly_files.keys():
            for r in input_reference_assembly_files.keys():
                smk_targets = smk_targets + [ "comparisons/" + r + "/" + s + ".pdf" ]

    if workflow_to_run == "compare-mummer":
        workflow_to_run = "compare-syri"
        for s in input_sample_assembly_files.keys():
            for r in input_reference_assembly_files.keys():
                smk_targets = smk_targets + [ "02_mummer_results/" + r + "/" + s + ".coords" ]

    if workflow_to_run == "evaluate-nanopore-reads":
        smk_targets = smk_targets + [ "nanopore_read_stats/{}".format(key) for key in input_nanopore_files ]

    if workflow_to_run == "evaluate-inspector":
        smk_targets = smk_targets + [ "inspector_assembly_evaluation/{}".format(key) for key in input_assembly_files ]

    if workflow_to_run == "evaluate-coverage":
        smk_targets = smk_targets + [ "evaluate/coverage_plots/nanopore/{}".format(key) for key in input_assembly_files ]

    if workflow_to_run == "evaluate-breseq-coverage":
        smk_targets = smk_targets + [ "evaluate/coverage_plots/breseq_nanopore/{}".format(key) for key in input_assembly_files ]

    if workflow_to_run == "evaluate-isescan":
        smk_targets = smk_targets + [ "evaluate/isescan/{}.csv".format(key) for key in input_assembly_files ]
    
    if workflow_to_run == "evaluate-soft-clipping":
        smk_targets = smk_targets + [ "evaluate/soft_clipping_summary/nanopore/{}_soft_clipping_summary.csv".format(key) for key in input_assembly_files ]
        config_options_list.append("brefito_package_path=" + str(brefito_package_path))

    if workflow_to_run == "evaluate-redotable":
        smk_targets = smk_targets + [ "evaluate/dot_plot/{}.svg".format(key) for key in input_assembly_files ]
        #config_options_list.append("brefito_package_path=" + str(brefito_package_path))


    #################################################
    ### trycycler trifecta
    #################################################

    if workflow_to_run == "trycycler-assemble":
        smk_targets = [ "05_trycycler/" + d + "/done" for d in input_nanopore_files.keys() ]
        #resource_options_list = resource_options_list + ["necats=4"]

    if workflow_to_run == "trycycler-reconcile":
        input_files=glob.glob("05_trycycler/*/cluster_*")
        for this_input_file in input_files:
            smk_targets.append(os.path.join(this_input_file, "2_all_seqs.fasta"))

    if workflow_to_run == "trycycler-consensus":
        smk_targets = [ "assemblies/" + d + ".fasta" for d in input_nanopore_files.keys() ]

    #################################################
        
    ### <---- Commands that haven't yet been updated above

    smk_file_path = os.path.join(rules_path, workflow_to_run + ".smk")
    target_options = ["-s", smk_file_path]

    config_options_list.append("references=" + references_argument)

    if samples_to_run != None and len(samples_to_run)>0:
        config_options_list.append("samples=" + "_,_".join(samples_to_run))


    config_options = []
    if len(config_options_list) > 0:
        config_options =  ["--config"]
        for i in config_options_list:
            si = i.split('=')
            config_options = config_options + [si[0] + '="' + si[1] + '"']

    if workflow_to_run=="search-blast":
        config_options = config_options + ["QUERY_FILE_PATH" + '="' + search_query_option + '"']

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
    print("RUNNING SNAKEMAKE COMMAND")
    print()
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

    if workflow_to_run == "polish-breseq":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "breseq")
    if workflow_to_run == "polish-polypolish":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "polypolish")
    elif workflow_to_run == "polish-polca":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "polca")
    elif workflow_to_run == "polish-medaka":
        copy_and_rename_assemblies(input_assembly_files.values(), "polished", "medaka")
    elif workflow_to_run == "normalize-assemblies":
        copy_and_rename_assemblies(input_assembly_files.values(), "normalized", "normalized")

def check_command_with_references(workflow_to_run, test_command_prefix):
    return_dict = { 'matched' : False }

    if workflow_to_run == test_command_prefix:
        return_dict['matched'] = True
        return_dict['references'] = None
        return_dict['workflow_to_run'] = test_command_prefix
    
    if workflow_to_run.startswith(test_command_prefix + "-"):
        return_dict['matched'] = True
        return_dict['references'] = workflow_to_run[len(test_command_prefix + "-"):]
        return_dict['workflow_to_run'] = test_command_prefix

    return (return_dict)

def check_command_list_with_references(workflow_to_run, test_command_prefix_list):

    for p in test_command_prefix_list:
        this_return_dict = check_command_with_references(workflow_to_run, p)
        if (this_return_dict['matched']):
            return (this_return_dict)

    return ( { 'matched' : False } )


if __name__ == "__main__":
    main()
