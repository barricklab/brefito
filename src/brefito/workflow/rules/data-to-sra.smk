try: sample_info
except NameError: 
    include: "load-sample-info.smk"

rule data_to_sra:
    output:
        "sra.csv"
    threads: 1
    run:
        with open(output[0], 'w') as f:
            for s in sorted(sample_info.get_sample_list()):
                f.write(",".join( [s] + sample_info.get_illumina_PE_read_list(s) + sample_info.get_illumina_SE_read_list(s) + sample_info.get_nanopore_read_list(s)) + "\n")