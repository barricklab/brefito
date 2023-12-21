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

print("HERHEHREE")

def all_local_paths():
    return(sample_info.get_local_path_list())

rule download_all:
    input: all_local_paths()

rule download_file:
    output:
        output_path = "{download_type}/{sample}",
    params:
        URL=lambda wildcards: sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample),
        bookmark = BOOKMARK,
        download_path = "{download_type}/temp_{sample}",
        lftp_commands_file = "{download_type}/temp_{sample}.lftp.commands"
    threads: 1
    wildcard_constraints:
        download_type="(references|illumina_reads|nanopore_reads)"
    resources:
        # This is an invented resource to prevent opening too many download connections at once!
        connections=1
    shell:
        """
        rm -f {params.download_path}
        echo 'cd "{REMOTE_BASE_PATH}"' > {params.lftp_commands_file}
        echo 'get "{params.URL}" -o "{params.download_path}"' >> {params.lftp_commands_file}
        #cat {params.lftp_commands_file}
        lftp {params.bookmark} < {params.lftp_commands_file}
        rm {params.lftp_commands_file}
        mv {params.download_path} {output.output_path}
        """