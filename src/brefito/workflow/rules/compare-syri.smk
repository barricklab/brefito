rule sort_reference_assemblies:
    input:
        "references/{reference}.fasta",
    output:
        "01_normalized_assemblies/references/{reference}.fasta"
    shell:
        "normalize_assembly -r {input} -i {input} -o {output} -s"


rule normalize_assemblies:
    input:
        reference = "01_normalized_assemblies/references/{reference}.fasta",
        sample = "samples/{sample}.fasta"
    output:
        "01_normalized_assemblies/samples/{reference}/{sample}.fasta"
    shell:
        "normalize_assembly -r {input.reference} -i {input.sample} -o {output} -s -c -x"

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
        mkdir -p 02_mummer_results/{wildcards.reference}
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