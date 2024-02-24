#! /usr/bin/env python3

# This script:
#   Renames contigs to match or be in order 
#   Sorts contigs in reverse order by length
#   Reindexes (shifts) the sequence in the input file to a requested subsequence

# ---- Imports ----
import argparse
import os
import importlib.util
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def iterate_files(directory):
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        if os.path.isfile(filepath):
            yield filepath

def main():
    # Get arguments
    parser = argparse.ArgumentParser(description="Extract ")
    parser.add_argument("-i", "--input", help="Input sequence file", action="store", required=True)
    parser.add_argument("-o", "--output", help="Output sequence file", action="store", default=None)
    parser.add_argument("-t", "--file-type", help="Input sequence file type [fasta]", action="store", default="fasta")


    args = parser.parse_args()

    mobile_element_seq_dict = {}

    # Load reference file
    
    reference_seqs = []

    file_path = args.input
    print("Reference File: ", file_path)
    for record in SeqIO.parse(file_path, args.file_type):
        num = 0
        for feature in record.features:
            if feature.type == "mobile_element":
                mobile_element_name = feature.qualifiers["name"][0]
                mobile_element_seq = feature.location.extract(record)
                num = num + 1
                mobile_element_seq.id = record.id + "_" + str(num)
                mobile_element_seq.name = record.id + "_" + str(num)
                #print(mobile_element_seq)

                mobile_element_seq.annotations['molecule_type'] = record.annotations['molecule_type']
                mobile_element_seq_dict[mobile_element_seq.seq] = mobile_element_seq


    records_to_write = []

    for seq in mobile_element_seq_dict.keys():

        records_to_write.append( mobile_element_seq_dict[seq] )

    SeqIO.write(records_to_write, args.output, "genbank")



if __name__ == '__main__':
    main()

