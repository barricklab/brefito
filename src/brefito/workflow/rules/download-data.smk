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

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

def remote_path_to_shell_command(remote_path, download_path):
    #print("Downloading: " + remote_path)

    split_remote_path = remote_path.split('://', 1)
    method = None
    bookmark = ''
    if (len(split_remote_path) == 1):
        sys.exit ('  FAILED: URL with :// not found.') 
    else: 
        # lftp@utbox://path/from/bookmark
        split_protocol = split_remote_path[0].split('@', 1)
        if (split_protocol[0] == 'lftp'):
            method = 'lftp'
            if (len(split_protocol) == 2): 
                bookmark = split_protocol[1]
            remote_path = split_remote_path[1]
        # ftp://server.com/path/to/file
        else:
            method = 'wget'

    if method == None:
        sys.exit ('  FAILED: Could not determine download type.')
    
    shell_command = ''
    if method == 'wget':
        shell_command = "wget -O \"{}\" {}".format(download_path, remote_path)
    elif method == 'lftp':
        #echo 'cd "{REMOTE_BASE_PATH}"' > {params.lftp_commands_file}
        #echo 'get "{params.URL}" -o "{params.download_path}"' >> {params.lftp_commands_file}
        shell_command = "echo 'get \"{}\" -o \"{}\"' | lftp {} ".format(remote_path, download_path, bookmark)

    #print(shell_command)
    return(shell_command)

def all_local_paths():
    return(sample_info.get_local_path_list())

rule download_all:
    input: all_local_paths()

rule download_file:
    output:
        output_path = "{download_type}/{sample}",
    params:
        remote_path = lambda wildcards: sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample),
        download_path = "{download_type}/temp_{sample}",
        shell_command = lambda wildcards:  remote_path_to_shell_command(sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample), wildcards.download_type + "/temp_" + wildcards.sample)
    threads: 1
    wildcard_constraints:
        download_type="(references|illumina_reads|nanopore_reads)"
    resources:
        # This is an invented resource to prevent opening too many download connections at once!
        connections=1
    conda:
        "../envs/download.yml"
    shell:
        """
        # Clear any stale output
        rm -f {params.download_path}
        {params.shell_command}
        mv {params.download_path} {output.output_path}
        """