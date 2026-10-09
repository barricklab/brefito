try: sample_info
except NameError:
    include: "load-sample-info.smk"

rule all_evaluate_isescan:
    input:
        ["evaluate/isescan/" + s + ".csv" for s in sample_info.get_sample_list()]
    default_target: True

rule evaluate_isescan:
    input:
        "assemblies/{sample}.fasta",
    output:
        dir = temp(directory("intermediates/isescan/{sample}")),
        file = "evaluate/isescan/{sample}.csv" 
    log:
        "logs/{sample}/isescan.log"
    conda:
        "../envs/isescan.yml"
    threads: 8
    shell:
        """
        isescan.py --nthread {threads} --seqfile {input} --output {output.dir} > {log} 2>&1
        mv {output.dir}/assemblies/{wildcards.sample}.fasta.csv {output.file}
        """