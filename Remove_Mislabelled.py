#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
from tqdm import tqdm

"""
Remove mislabelled sequences identified by Sativa algorithm.

All the sequences identified by Sativa are removed from
the genbank file.
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
        description="Remove mislabelled sequences.")
    parser.add_argument("gb", type=argparse.FileType("r"),
                   help="genbank file (.gb)")
    parser.add_argument("mis", type=argparse.FileType("r"),
                    help="file with mislabelled sequences" \
                    " resulted after running sativa (.mis)")
    
    return parser.parse_args()

def parse_mis(mis):
    """Read file.mis and parse all sequences which need
    to be removed.

    This function parses the acession numbers in file.mis
    that should be removed from the genbank file.

    Parameters
    -------------
    mis : io.TextIOWrapper
        Sativa mislabelled sequences file

    Returns
    -------------
    acession_numbers: list
        List of acession numbers
    """

    acession_numbers = {m.split("\t")[0]:m.split("\t")[4] for m in mis
                        if not m.startswith(";")}

    return acession_numbers

def remove_mislabelled(gb, acession_numbers):
    """Remove mislabelled sequences from genbank file

    This function removes from the genbank file all
    sequences found by Sativa algorithm.

    Parameters
    -------------
    gb : io.TextIOWrapper
        Genbank file

    acession_numbers : list
        List of acession numbers to be removed

    """

    gb_all_data = SeqIO.index(gb.name, "genbank")

    base = os.path.basename(gb.name)
    name = os.path.splitext(base)[0]

    gb_m_name = "{}_m.gb".format(name)
    with open(gb_m_name, "w") as gb_c:
        for acession in tqdm(gb_all_data, desc="Removing mislabelled from genbank"):
            if not acession in acession_numbers:
                print(gb_all_data.get_raw(acession).decode(),
                            file=gb_c, end="")


if __name__ == "__main__":
    args = getArguments()
    acession_numbers = parse_mis(args.mis)
    remove_mislabelled(args.gb, acession_numbers)