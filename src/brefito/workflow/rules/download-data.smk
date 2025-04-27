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

import os
import uuid

# Use multiple threads for SRA downloads to speed up gzipping
def get_num_threads(remote_path):
    if (remote_path.startswith("sra://")):
        return 6
    return 1

def get_num_connections(remote_path):
    if (remote_path.startswith("sra://")):
        return 0
    return 1

def remote_path_to_shell_command(remote_path, download_path, output_path):
    #print("Downloading: " + remote_path)

    #Save b/c it can be modified below and we need to keep the prefix for some functionality
    original_remote_path = remote_path

    split_remote_path = remote_path.split('://', 1)
    method = None
    bookmark = ''
    if (len(split_remote_path) == 1):
        #Check that the file exists to give better error messages
        if not os.path.exists(remote_path):
            sys.exit ('  FAILED: Local file does not exist: ' + remote_path)

        #create a symbolic link...
        shell_command = "echo \"Creating local symbolic link for file that does not have a download URL: " + remote_path + "\""

        if os.path.isabs(remote_path):
            shell_command = shell_command + " && ln -s " + remote_path + " " + output_path
        else: #Add ../ to move down a level if it is a local path
            shell_command = shell_command + " && ln -s ../" + remote_path + " " + output_path

        #print(shell_command)
        return(shell_command)
    else:
        # lftp@utbox://path/from/bookmark
        split_protocol = split_remote_path[0].split('@', 1)
        if (split_protocol[0] == 'lftp'):
            method = 'lftp'
            if (len(split_protocol) == 2):
                bookmark = split_protocol[1]
            remote_path = split_remote_path[1]
        # ncbi://accession
        elif (split_protocol[0] == 'ncbi'):
            method = 'ncbi'
            remote_path = split_remote_path[1]
        # ftp://server.com/path/to/file or https://server.com/path/to/file etc.

        # sra://SRRXXXX sra accession
        elif (split_protocol[0] == 'sra'):
            method = 'sratools'
            remote_path = split_remote_path[1]

        else:
            method = 'wget'

    if method == None:
        sys.exit ('  FAILED: Could not determine download type.')

    shell_command = ''

    if method == 'wget':
        shell_command = "wget -O \"{}\" {}".format(download_path, remote_path)
        ## Move the temp download path to the final path
        shell_command = shell_command + " && mv " + download_path + " " + output_path

    elif method == 'lftp':
        #echo 'cd "{REMOTE_BASE_PATH}"' > {params.lftp_commands_file}
        #echo 'get "{params.URL}" -o "{params.download_path}"' >> {params.lftp_commands_file}
        shell_command = "echo 'get \"{}\" -o \"{}\"' | lftp {} ".format(remote_path, download_path, bookmark)
        ## Move the temp download path to the final path
        shell_command = shell_command + " && mv " + download_path + " " + output_path

    elif method == 'ncbi':
        shell_command = "sleep 1; esearch -db nucleotide -query {} | efetch -format genbank > {}".format(remote_path, download_path)
        ## Move the temp download path to the final path
        shell_command = shell_command + " && mv " + download_path + " " + output_path

    elif method == 'sratools':

        ## Downloading is handled by other rules for efficiency/speed. Here we only move them to final locations

        if output_path.endswith("SE.fastq.gz"): #this is SE reads
            shell_command = f"mv sra-downloads/{remote_path}.SE.fastq.gz {output_path}"
        elif output_path.endswith("R1.fastq.gz"): #this is R1 of PE data
            shell_command = f"mv sra-downloads/{remote_path}.R1.fastq.gz {output_path}"
        elif output_path.endswith("R2.fastq.gz"): #this is R2 of PE data
            shell_command = f"mv sra-downloads/{remote_path}.R2.fastq.gz {output_path}"
        elif "nanopore" in output_path: #this is nanopore data
            shell_command = f"mv sra-downloads/{remote_path}.fastq.gz {output_path}"
        return shell_command

    #print(shell_command)
    return(shell_command)

def all_local_paths():
    return(sample_info.get_local_path_list())

rule download_all:
    input: all_local_paths()

# We need to download the SRA stuff for some inputs
def get_sra_path_from_local_path(local_path, file_name):
    remote_path = sample_info.get_remote_path_from_local_path(local_path)
    if not remote_path.startswith("sra://"):
        return []

    #print(f"sra-downloads/{file_name}")
    return(f"sra-downloads/{file_name}")

rule download_file:
    output:
        output_path = "{download_type}/{sample}",
    input:
        input_sra_path = lambda wildcards: get_sra_path_from_local_path(wildcards.download_type + "/" + wildcards.sample, wildcards.sample),
    params:
        remote_path = lambda wildcards: sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample),
        download_path = "{download_type}/temp_{sample}",
        shell_command = lambda wildcards: remote_path_to_shell_command(sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample), wildcards.download_type + "/temp_" + wildcards.sample, wildcards.download_type + "/" + wildcards.sample)
    threads: 1
    wildcard_constraints:
        download_type="(references|illumina-reads|nanopore-reads)"
    resources:
        # This is an invented resource to prevent opening too many download connections at once!
        # We calculate this so we can use 0 here for SRA downloads handled by a separate function
        connections=lambda wildcards: get_num_connections(sample_info.get_remote_path_from_local_path(wildcards.download_type + "/" + wildcards.sample))
    conda:
        "../envs/download.yml"
    shell:
        """
        # Clear any stale output
        rm -f {params.download_path}
        {params.shell_command}
        """

rule download_sra:
    output:
        temp("sra-downloads/{run_accession}/{run_accession}.sra")
    threads: 1
    resources:
        connections=1
    conda:
        "../envs/download.yml"
    shell:
        """
        mkdir -p sra-downloads
        cd sra-downloads
        prefetch {wildcards.run_accession}
        """

ruleorder: dump_paired_sra > dump_single_end_sra > dump_nanopore_sra

rule dump_paired_sra:
    input:
        "sra-downloads/{run_accession}/{run_accession}.sra"
    output:
        output_R1_path = temp("sra-downloads/{run_accession}.R1.fastq"),
        output_R2_path = temp("sra-downloads/{run_accession}.R2.fastq")
    threads: 6
    conda:
        "../envs/download.yml"
    shell:
        """
        cd sra-downloads
        fasterq-dump -e {threads} {wildcards.run_accession}
        mv {wildcards.run_accession}_1.fastq {wildcards.run_accession}.R1.fastq
        mv {wildcards.run_accession}_2.fastq {wildcards.run_accession}.R2.fastq
        """

rule dump_single_end_sra:
    input:
        "sra-downloads/{run_accession}/{run_accession}.sra"
    output:
        temp("sra-downloads/{run_accession}.SE.fastq")
    threads: 6
    conda:
        "../envs/download.yml"
    conda:
        "../envs/download.yml"
    shell:
        """
        cd sra-downloads
        fasterq-dump -e {threads} {wildcards.run_accession}
        mv {wildcards.run_accession}.fastq {wildcards.run_accession}.SE.fastq
        """

rule dump_nanopore_sra:
    input:
        "sra-downloads/{run_accession}/{run_accession}.sra"
    output:
        temp("sra-downloads/{run_accession}.fastq")
    threads: 6
    conda:
        "../envs/download.yml"
    conda:
        "../envs/download.yml"
    shell:
        """
        cd sra-downloads
        fasterq-dump -e {threads} {wildcards.run_accession}
        """

rule gzip_sra:
    input:
        "sra-downloads/{file_name}.fastq"
    output:
        temp("sra-downloads/{file_name}.fastq.gz")
    threads: 6
    conda:
        "../envs/download.yml"
    shell:
        """
        pigz -f -p {threads} {input}
        """