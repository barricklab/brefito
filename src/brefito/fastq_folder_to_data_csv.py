#!/usr/bin/env python3

def main():

    import glob
    import os.path
    import argparse
    import re
    import csv

    # What command did we choose
    parser = argparse.ArgumentParser(
                        prog='fastq_folder_to_data_csv_folder',
                        description='automatically create a data.csv file for brefito from an input folder of FASTQ files',
                        epilog='')

    parser.add_argument('-i', '--input', default='.', type=str, help='Path to input folder containing FASTQ files')
    parser.add_argument('-o', '--output', default='data.csv', type=str, help='Path to output data.csv')
    parser.add_argument('-p', '--prefix', default='', type=str, help='Prefix to add to paths (remote location)')

    args = parser.parse_args()

    input_path = args.input
    output_path = args.output
    output_prefix = args.prefix

    
    print("Looking for *.fastq.gz files in: " + input_path)

    in_file_ending = "fastq.gz"

    fastq_files=glob.glob(os.path.join(input_path, "*."+in_file_ending))
    fastq_files=sorted(fastq_files)
    input_files = []
    for this_input_file in fastq_files:
        this_file_name=os.path.basename(this_input_file)
        #print(this_file_name)
        #short_name = re.findall(r'(.+)\.' + re.escape(in_file_ending), this_file_name)
        input_files += [this_file_name]

    #print("\n".join(input_files))

    # Group / fix paired reads that end in _1/2
    line_list = [ ['sample', 'type', 'setting'] ]
    i=0
    while i<len(input_files):

        written = False
        if i!=len(input_files)-1:
            j = i+1

            short_name_i = input_files[i].replace('_1.fastq.gz', '')
            short_name_j = input_files[j].replace('_2.fastq.gz', '')

            if (short_name_i == short_name_j):
                line_list.append([short_name_i, 'illumina-PE', output_prefix +  short_name_i + '{1|2}.fastq.gz'])
                i+=1
                written = True

        if not written:
            short_name = input_files[i].replace('.fastq.gz', '')
            line_list.append([short_name, 'illumina-SE', output_prefix + input_files[i]])

        i += 1

    #print(line_list)

    if len(line_list)==1:
        print("No input files found. No data.csv written.")
        return

    print("Writing data.csv file: " + output_path)
    print("     with path prefix: " + output_prefix)

    with open(output_path, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)
        writer.writerows(line_list)

if __name__ == "__main__":
    main()