import glob
import os.path

INPUT_SAMPLE_FILES=glob.glob("input/samples/*.fasta")

INPUT_REFERENCE_FILES=glob.glob("input/references/*.fasta")
print(INPUT_REFERENCE_FILES)
OUTPUT_PDF_FILES = []
SYRI_OUTPUT_FILES = []
for this_reference_file in INPUT_REFERENCE_FILES:
    this_file_name=os.path.basename(this_reference_file)
    reference_name = re.findall("(.+)\.fasta", this_file_name)


    for this_sample_file in INPUT_SAMPLE_FILES:
        this_sample_file_name=os.path.basename(this_sample_file)
        sample_name = re.findall("(.+)\.fasta", this_sample_file_name)
        SYRI_OUTPUT_FILES.append(os.path.join("02_syri", reference_name[0], sample_name[0],  "syri.out"))

        OUTPUT_PDF_FILES.append(os.path.join("output", reference_name[0], sample_name[0] + ".pdf"))

print(OUTPUT_PDF_FILES)

print(SYRI_OUTPUT_FILES)



    

rule all:
    input: OUTPUT_PDF_FILES

rule compare_mummer:
    input:
        reference = "input/references/{reference}.fasta",
        sample = "input/samples/{sample}.fasta"
    output:
        delta = "01_mummer/{reference}/{sample}.delta",
        filtered_delta = "01_mummer/{reference}/{sample}.filtered.delta",
        coords = "01_mummer/{reference}/{sample}.coords"
    conda:
        "../envs/mummer_syri_plotsr.yml"
    shell:
        """
        nucmer --maxmatch -c 100 -b 500 -p 01_mummer/{wildcards.reference}/{wildcards.sample} {input.reference} {input.sample}
        delta-filter -m -i 90 -l 100 {output.delta} > {output.filtered_delta}
        show-coords -THrd {output.filtered_delta} > {output.coords}
        """

rule create_genomes_tsv:
        input:
            reference_fasta_file = "input/references/{reference}.fasta",
            sample_fasta_files = "input/samples/{sample}.fasta"
        output:
            "03_plotsr/{reference}/{sample}/genomes.tsv"
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
        reference = "input/references/{reference}.fasta",
        sample = "input/samples/{sample}.fasta",
        coords = "01_mummer/{reference}/{sample}.coords",
        filtered_delta = "01_mummer/{reference}/{sample}.filtered.delta"
    output:
        "02_syri/{reference}/{sample}/syri.out"
    conda:
        "../envs/mummer_syri_plotsr.yml"
    shell:
        """
        mkdir -p 02_syri/{wildcards.reference}/{wildcards.sample}
        syri --nosnp --dir 02_syri/{wildcards.reference}/{wildcards.sample} -c {input.coords} -d {input.filtered_delta} -r {input.reference} -q {input.sample}
        """

def find_available_files(wildcards):
   from glob import glob
   path = "02_syri/{reference}/{sample}/syri.out"
   files = glob(path.format(reference=wildcards.reference, sample="*"))
   return files

rule run_plotsr:
    input:
        syri_out = "02_syri/{reference}/{sample}/syri.out",
        genomes_tsv = "03_plotsr/{reference}/{sample}/genomes.tsv"
    output:
        "output/{reference}/{sample}.pdf"
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
        plotsr $srarg --genomes {input.genomes_tsv} -H 4 -W 8 -s 100 -o {output}
        """