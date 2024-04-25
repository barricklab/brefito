# Generates compare tables from all sets of samples that share references
include: "predict-mutations-breseq.smk"

try: sample_info
except NameError: 
    include: "load-sample-info.smk"

# Divide up the samples into sets that share the exact same reference files.

# Hashed by a key that concatenates the references using '|||''...
sample_lists_by_references = {}

for s in sample_info.get_sample_list():
    concatenated_reference_key = "|||".join(sample_info.get_reference_list(s))
    if not concatenated_reference_key in sample_lists_by_references:
        sample_lists_by_references[concatenated_reference_key] = []
    sample_lists_by_references[concatenated_reference_key].append(s)
#print(sample_lists_by_references)

num_sample_sets = len(sample_lists_by_references)
sorted_keys = sorted(sample_lists_by_references, key=lambda x: len(sample_lists_by_references[x]))

def output_name_from_key_index(key_index):
    if (len(sorted_keys)==0):
        return 'all'
    return str(key_index)

def output_name_to_key_index(output_name):
    shorten_output_name = re.sub("^compare_*", "", output_name)
    shorten_output_name = re.sub("\.html$", "", shorten_output_name)
    #should be a number or all at this point
    if shorten_output_name == 'all':
        return 0
    return int(shorten_output_name)-1

def get_gd_list_from_output_name(output_name):
    key_index = output_name_to_key_index(output_name)
    sample_list = sample_lists_by_references[sorted_keys[key_index]]
    gd_list = ["breseq-" + sample_info.get_reference_prefix() + "/gd/" + s + ".gd" for s in sample_list]
    return(gd_list)

def get_reference_files_list_from_output_name(output_name):
    key_index = output_name_to_key_index(output_name)
    reference_key = sorted_keys[key_index]
    reference_list = reference_key.split('|||')
    return(reference_list)


rule all_compare_mutations_breseq:
    input:
        ["breseq-" + sample_info.get_reference_prefix() + "/compare_" + output_name_from_key_index(i+1)  + ".html" for i in range(num_sample_sets)]
    default_target: True

rule compare_mutations_breseq:
    input:
        gd_files = lambda wildcards: get_gd_list_from_output_name(wildcards.output_name),
        reference_files = lambda wildcards: get_reference_files_list_from_output_name(wildcards.output_name)
    output:
        "breseq-" + sample_info.get_reference_prefix() + "/compare_{output_name}.html"
    log: 
        "logs/breseq-compare-" + sample_info.get_reference_prefix() + "-{output_name}.log"
    conda:
        "../envs/breseq.yml"
    params:
        reference_arguments = lambda wildcards: ["-r " + r for r in get_reference_files_list_from_output_name(wildcards.output_name)]
    threads: 1
    shell:
        """
        gdtools COMPARE -o {output} {params.reference_arguments} {input.gd_files} > {log} 2>&1
        """