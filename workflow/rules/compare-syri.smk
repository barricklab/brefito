# import glob
# import os.path

# INPUT_REFERENCE_FILES=glob.glob("input/references/*.fasta")
# print(INPUT_REFERENCE_FILES)

# INPUT_SAMPLE_FILES=glob.glob("input/samples/*.fasta")
# print(INPUT_SAMPLE_FILES)

# OUTPUT_PDF_FILES = []
# SYRI_OUTPUT_FILES = []
# for this_reference_file in INPUT_REFERENCE_FILES:
#     this_file_name=os.path.basename(this_reference_file)
#     reference_name = re.findall("(.+)\.fasta", this_file_name)


#     for this_sample_file in INPUT_SAMPLE_FILES:
#         this_sample_file_name=os.path.basename(this_sample_file)
#         sample_name = re.findall("(.+)\.fasta", this_sample_file_name)
#         SYRI_OUTPUT_FILES.append(os.path.join("02_syri", reference_name[0], sample_name[0],  "syri.out"))

#         OUTPUT_PDF_FILES.append(os.path.join("output", reference_name[0], sample_name[0] + ".pdf"))

# print(OUTPUT_PDF_FILES)


NORMALIZE_ASSEMBLIES_SCRIPT_PATH = os.path.join(workflow.current_basedir, "..", "scripts", "normalize_assembly.py")

# rule all:
#     input: OUTPUT_PDF_FILES

rule sort_reference_assemblies:
    input:
        "references/{reference}.fasta",
    output:
        "01_normalized_assemblies/references/{reference}.fasta"
    conda:
        "../envs/biopython.yml"
    shell:
        "{NORMALIZE_ASSEMBLIES_SCRIPT_PATH} -r {input} -i {input} -o {output} -s"


rule normalize_assemblies:
    input:
        reference = "01_normalized_assemblies/references/{reference}.fasta",
        sample = "samples/{sample}.fasta"
    output:
        "01_normalized_assemblies/samples/{reference}/{sample}.fasta"
    conda:
        "../envs/biopython.yml"
    shell:
        "{NORMALIZE_ASSEMBLIES_SCRIPT_PATH} -r {input.reference} -i {input.sample} -o {output} -s -c -x"

rule compare_mummer:
    input:
        reference = "01_normalized_assemblies/references/{reference}.fasta",
        sample = "01_normalized_assemblies/samples/{reference}/{sample}.fasta"
    output:
        delta = "02_mummer_results/{reference}/{sample}.delta",
        filtered_delta = "02_mummer_results/{reference}/{sample}.filtered.delta",
        coords = "02_mummer_results/{reference}/{sample}.coords",
    conda:
        "../envs/mummer_syri_plotsr.yml"
    shell:
        """
        nucmer --maxmatch -c 100 -b 500 -p 02_mummer_results/{wildcards.reference}/{wildcards.sample} {input.reference} {input.sample}
        delta-filter -m -i 90 -l 100 {output.delta} > {output.filtered_delta}
        show-coords -THrd {output.filtered_delta} > {output.coords}
        """

rule create_genomes_tsv:
        input:
            reference_fasta_file = "01_normalized_assemblies/references/{reference}.fasta",
            sample_fasta_files = "01_normalized_assemblies/samples/{reference}/{sample}.fasta"
        output:
            "04_plotsr_results/{reference}/{sample}/genomes.tsv"
        shell:
            """
            for i in {input.reference_fasta_file} {input.sample_fasta_files} 
            do
                #echo $i
                [[ $i =~ .+/(.+).fasta$ ]] && echo -e "$i\t${{BASH_REMATCH[1]}}\tft:fa;lw:1" >> {output}
            done
            """

rule run_syri:
    input:
        reference = "01_normalized_assemblies/references/{reference}.fasta",
        sample = "01_normalized_assemblies/samples/{reference}/{sample}.fasta",
        coords = "02_mummer_results/{reference}/{sample}.coords",
        filtered_delta = "02_mummer_results/{reference}/{sample}.filtered.delta"
    output:
        syri_out = "03_syri_results/{reference}/{sample}/syri.out",
        syri_dir = directory("03_syri_results/{reference}/{sample}")
    conda:
        "../envs/mummer_syri_plotsr.yml"
    shell:
        """
        mkdir -p {output.syri_dir}
        syri --nosnp --dir {output.syri_dir} -c {input.coords} -d {input.filtered_delta} -r {input.reference} -q {input.sample}
        """

def find_available_files(wildcards):
   from glob import glob
   path = "03_syri_results/{reference}/{sample}/syri.out"
   files = glob(path.format(reference=wildcards.reference, sample="*"))
   return files

rule run_plotsr:
    input:
        syri_out = "03_syri_results/{reference}/{sample}/syri.out",
        genomes_tsv = "04_plotsr_results/{reference}/{sample}/genomes.tsv"
    output:
        "comparisons/{reference}/{sample}.pdf"
    log: 
        "logs/plotsr_{reference}_{sample}.log"
    conda:
        "../envs/mummer_syri_plotsr.yml"
    shell:
        """
        srarg=""
        for i in {input.syri_out} 
        do
            echo $i
            srarg="$srarg --sr $i"
        done
        echo $srarg
        plotsr $srarg --genomes {input.genomes_tsv} -H 4 -W 8 -s 100 -o {output} --lf {log}
        """