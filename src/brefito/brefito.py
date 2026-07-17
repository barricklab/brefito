#!/usr/bin/env python3

import glob
import os
import os.path
import argparse
import re
import shutil
import subprocess
import hashlib
import json
import sys

def main():

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

    parser.add_argument('-p', '--path', default='.', type=str,
                         help='Working directory in which to run the workflow. brefito '
                              'changes to this directory before discovering input files '
                              'and running Snakemake, so all sample/reference/output '
                              'paths are resolved relative to it. Default: current directory')

    #An additional way to specify these
    parser.add_argument('-r', '--references', type=str, help='Path to use for reference files. Default: references')
    parser.add_argument('-d', '--data', default=None,
                        help='path to CSV with sample and file information (default: data.csv). '
                             'If specified explicitly, the file must exist.')

    #Snakemake passthroughs
    parser.add_argument('--cores', type=int, default=0, help='--cores argument passed through to Snakemake (0 = all)') # 0 means "all"
    parser.add_argument('--config', action='append', default=[], help='--config argument passed through to Snakemake. Individual workflows support different settings.')
    parser.add_argument('--resources', action='append', default=[], help='--resources argument passed through to Snakemake. Individual workflows support different settings.')
    parser.add_argument('--rerun-incomplete', action='store_true', help='argument passed through to Snakemake') 
    parser.add_argument('--rerun-triggers', nargs='+', default=None, help='argument passed through to Snakemake. Takes one or more space-separated triggers under a single flag (e.g. "--rerun-triggers mtime code input"), matching Snakemake\'s own syntax. Providing it replaces the default set (mtime, code, input) rather than adding to it. Add "params" to also re-run when a rule\'s params change (e.g. to force re-download of references).')
    parser.add_argument('--unlock', action='store_true', help='argument passed through to Snakemake') 
    parser.add_argument('--pick-lock', action='store_true', help='run Snakemake --unlock and then immediately run Snakemake')
    parser.add_argument('--reinstall', action='append', default=[], metavar='ENV',
                         help='Delete the installed Snakemake conda environment matching '
                              'workflow/envs/ENV or workflow/envs/ENV.yml, forcing Snakemake '
                              'to reinstall it on this run. Can be specified multiple times. '
                              'Example: --reinstall breseq')
    parser.add_argument('--breseq-prerelease', action='store_true',
                         help='Use the prerelease breseq build (breseq-prerelease conda env) '
                              'instead of the stable release. Use --reinstall breseq-prerelease '
                              'to force-reinstall the prerelease env when a new build is pushed.')
    parser.add_argument('--notemp', action='store_true', help='argument passed through to Snakemake')
    parser.add_argument('--keep-going', action='store_true', help='argument passed through to Snakemake')
    parser.add_argument('--nolock', action='store_true', help='argument passed through to Snakemake')
    parser.add_argument('--forceall', action='store_true', help='argument passed through to Snakemake')
    parser.add_argument('-n', '--dry-run', action='store_true', help='argument passed through to Snakemake')

    # REQUIRED positional argument
    parser.add_argument('workflow', type=str)

     # OPTIONAL positional arguments
    parser.add_argument('samples', nargs='*', type=str)

    args = parser.parse_args()

    base_path = args.path
    # Run everything (input-file discovery, Snakemake, its .snakemake/conda dirs, and
    # output writing) in the requested directory by changing to it up front, so all
    # relative paths resolve consistently there rather than in the invocation CWD.
    if not os.path.isdir(base_path):
        print("Error: --path directory does not exist: " + base_path)
        sys.exit(1)
    os.chdir(base_path)
    workflow_to_run = args.workflow.lower()
    samples_to_run = args.samples

    config_options_list = args.config
    resource_options_list = args.resources
    references_argument = args.references

    # If the user explicitly specified a --data/-d file, it must exist. If they did
    # not, fall back to the default 'data.csv', which may legitimately be absent for
    # workflows that discover samples from the directory structure instead.
    if args.data is not None:
        if not os.path.isfile(args.data):
            print("Error: --data/-d file does not exist: " + args.data)
            sys.exit(1)
        data_file_name = args.data
    else:
        data_file_name = "data.csv"

    if args.pick_lock and args.unlock:
        print("Ignoring --unlock option because --pick-lock also provided.")

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


    snakemake_plus_common_options = ["snakemake", "--use-conda"]
    if args.cores == 0:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--cores", "all"]
    else:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--cores", str(args.cores)]

    rerun_triggers = args.rerun_triggers if args.rerun_triggers is not None else ['mtime', 'code', 'input']
    # Snakemake's --rerun-triggers uses nargs="+", so pass the whole list under a single
    # flag. (Emitting it once per trigger would make argparse keep only the last value,
    # silently dropping the other selected triggers such as 'code'.)
    snakemake_plus_common_options = snakemake_plus_common_options + ["--rerun-triggers"] + rerun_triggers

    if args.rerun_incomplete:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--rerun-incomplete"]
    if args.notemp:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--notemp"]
    if args.keep_going:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--keep-going"]
    if args.nolock:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--nolock"]
    if args.forceall:
        snakemake_plus_common_options = snakemake_plus_common_options + ["--forceall"]
    if args.dry_run:
        snakemake_plus_common_options = snakemake_plus_common_options + ["-n"]

    # What are appropriate targets for the workflow we are running?
    smk_targets = []

    if workflow_to_run in ["polish-breseq", "polish-polypolish", "polish-polca", "polish-medaka"]:
        print()
        print("Genome assembly files found (*.fasta) in " + assemblies_path)
        print()
        input_assembly_files = find_matching_input_files(assemblies_path, "fasta")
        for (k, v) in input_assembly_files.items(): print("    " + k + " : " + v)
        if (len(input_assembly_files.items()) == 0) : print("    " + "NONE FOUND")
        smk_targets = [ d + ".polished" for d in input_assembly_files.values() ]

    # Set this globally, putting it first means it can be overridden
    resource_options_list = ["connections=1"] + resource_options_list

    # Resurce so we can run just one of certain rules at a time
    resource_options_list = ["singleton=1"] + resource_options_list

    if workflow_to_run in ["check-soft-clipping", "curate-ltee-clones"]:
        config_options_list.append("brefito_package_path=" + str(brefito_package_path))

    # If not specified at command line or in workflow, set to default        
    if references_argument == None:
        references_argument = 'references'


    ## Commands that haven't yet been updated below --->

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

    config_options_list.append("data_csv=" + data_file_name)

    config_options_list.append("references=" + references_argument)

    if samples_to_run != None and len(samples_to_run)>0:
        config_options_list.append("samples=" + "_,_".join(samples_to_run))

    if args.breseq_prerelease:
        config_options_list.append("BRESEQ_PRERELEASE=True")

    config_options = []
    if len(config_options_list) > 0:
        config_options =  ["--config"]
        for i in config_options_list:
            # split once so a value that itself contains '=' is preserved, and so a
            # value-less flag (e.g. --config SKIP_PHISPY) becomes SKIP_PHISPY="1"
            # instead of crashing here on a missing si[1].
            si = i.split('=', 1)
            if len(si) == 1:
                config_options = config_options + [si[0] + '="1"']
            else:
                config_options = config_options + [si[0] + '="' + si[1] + '"']

    if workflow_to_run=="search-blast":
        config_options = config_options + ["QUERY_FILE_PATH" + '="' + search_query_option + '"']

    # This lets us replace defaults with user specified resource settings.
    # Snakemake resource names are lower case (e.g. 'connections'), so normalize
    # the key to lower case: this makes the passthrough case-insensitive and lets
    # a user-supplied 'CONNECTIONS=4' actually override the default 'connections=1'
    # instead of being sent as a separate, ignored resource.
    resource_options = []
    if len(resource_options_list) > 0:
        resource_options =  ["--resources"]
        resource_options_dict = {}
        for i in resource_options_list:
            si = i.split('=')
            resource_options_dict[si[0].lower()]=si[1]
        for k in resource_options_dict:
            resource_options = resource_options + [k + '=' + resource_options_dict[k]]

    command = snakemake_plus_common_options + target_options + smk_targets + config_options + resource_options

    if args.reinstall:
        envs_path = os.path.join(brefito_package_path, "workflow", "envs")
        reinstall_conda_envs(args.reinstall, envs_path, get_conda_envs_dir())

    print()
    print("RUNNING SNAKEMAKE COMMAND")
    print()

    snakemake_returncode = 0
    if args.unlock or args.pick_lock:
        print(" ".join(command + ["--unlock"]))
        subprocess.run(command + ["--unlock"])

    if not args.unlock:
        print(" ".join(command))
        snakemake_returncode = subprocess.run(command).returncode

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

    # Propagate Snakemake's exit status so failures (including a workflow that
    # requires sample information but found none) surface as a non-zero exit.
    sys.exit(snakemake_returncode)

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

def get_conda_envs_dir():
    conda_prefix = os.environ.get("SNAKEMAKE_CONDA_PREFIX")
    if conda_prefix:
        conda_prefix = os.path.expanduser(os.path.expandvars(conda_prefix))
        return os.path.abspath(conda_prefix)
    return os.path.join(os.path.abspath(".snakemake"), "conda")

def get_conda_platform():
    try:
        info = json.loads(subprocess.check_output(["conda", "info", "--json"], text=True))
        return info["platform"]
    except (subprocess.CalledProcessError, FileNotFoundError, KeyError, json.JSONDecodeError):
        print("Warning: could not determine conda platform (is conda on PATH?); ignoring any pin files.")
        return None

def candidate_env_hashes(yaml_file, envs_dir, platform):
    with open(yaml_file, "rb") as f:
        content = f.read()

    prefix = yaml_file
    if prefix.endswith(".yml") or prefix.endswith(".yaml"):
        prefix = prefix.rsplit(".", 1)[0]

    pin_content = None
    if platform:
        pin_file = prefix + f".{platform}.pin.txt"
        if os.path.isfile(pin_file):
            with open(pin_file, "rb") as f:
                pin_content = f.read()

    deploy_content = None
    deploy_file = prefix + ".post-deploy.sh"
    if os.path.isfile(deploy_file):
        with open(deploy_file, "rb") as f:
            deploy_content = f.read()

    # Mirrors Snakemake's Env._get_hash(include_location=True, include_container_img=True)
    current = hashlib.md5(usedforsecurity=False)
    current.update(os.path.realpath(envs_dir).encode())
    if deploy_content:
        current.update(deploy_content)
    if pin_content:
        current.update(pin_content)
    current.update(content)

    # Legacy fallbacks for envs built by older/different Snakemake versions.
    legacy_content_only = hashlib.md5(content).hexdigest()
    legacy_name_and_content = hashlib.md5(os.path.basename(yaml_file).encode() + content).hexdigest()

    return [current.hexdigest(), legacy_content_only, legacy_name_and_content]

def get_conda_executable():
    return "mamba" if shutil.which("mamba") else "conda"

def remove_installed_conda_env(yaml_file, envs_dir, platform):
    executable = get_conda_executable()
    for env_hash in candidate_env_hashes(yaml_file, envs_dir, platform):
        for candidate in (env_hash, env_hash + "_", env_hash[:8]):
            env_dir = os.path.join(envs_dir, candidate)
            if os.path.isdir(env_dir):
                subprocess.run([executable, "remove", "--all", "--yes", "--prefix", env_dir])
                for suffix in (".yaml", ".pin.txt", ".post-deploy.sh", ".env_setup_done"):
                    sidecar = env_dir + suffix
                    if os.path.isfile(sidecar):
                        os.remove(sidecar)
                return env_dir
    return None

def reinstall_conda_envs(reinstall_values, envs_path, envs_dir):
    if not reinstall_values:
        return

    platform = get_conda_platform()

    for value in reinstall_values:
        candidate_names = sorted({value, value + ".yml"})
        matches = [n for n in candidate_names if os.path.isfile(os.path.join(envs_path, n))]

        if not matches:
            print(f"--reinstall {value}: no environment file matching '{value}' or '{value}.yml' found in {envs_path}")
            continue

        for name in matches:
            yaml_file = os.path.join(envs_path, name)
            removed = remove_installed_conda_env(yaml_file, envs_dir, platform)
            if removed:
                print(f"--reinstall {value}: removed installed conda environment for {name} ({removed})")
            else:
                print(f"--reinstall {value}: {name} has no installed conda environment yet; it will be installed fresh")

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

if __name__ == "__main__":
    main()
