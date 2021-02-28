#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
from tqdm import tqdm

"""
Filter sequences by length and exclude all the outliers.

Find all sequences with length lower than threshold, and
remove them from the genbank file.
"""

def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    parser = argparse.ArgumentParser(
        description='Filter sequences by length and exclude all the outliers.')
    parser.add_argument('gb', type=argparse.FileType('r'),
                   help='genbank file .gb')
    parser.add_argument('-t', '--threshold', nargs='?', const=100, 
                    default=100, type = int,
                    help="Threshold for removing the sequences, default: 100")
    
    return parser.parse_args()

def filter_by_length(gb, threshold):
    """Filter sequences by length and exclude all the outliers

    This function removes all sequences with sequence length
    lower than the threshold.

    Parameters
    -------------
    gb : io.TextIOWrapper
        Genbank file

    threshold : int
        Threshold value
    """

    gb_all_data = SeqIO.index(gb.name, "genbank")
    gb_new_data = []

    base = os.path.basename(gb.name)
    name = os.path.splitext(base)[0]

    gb_l_name = "{}_l.gb".format(name)
    for acession in tqdm(gb_all_data, 
        desc="Removing sequences smaller than {}".format(threshold)):
        
        record = gb_all_data.get(acession)
        if len(record.seq) > threshold:
            gb_new_data.append(record)

    SeqIO.write(gb_new_data, gb_l_name, "genbank")


if __name__ == "__main__":
    args = getArguments()
    filter_by_length(args.gb, args.threshold)