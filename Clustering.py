#!/usr/bin/env python3
import os
import argparse
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from multiprocessing import cpu_count

"""Clustering sequences and removing redundancy from the file.

Removing redundancy by clustering sequences with vsearch.
""" 

def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    num_threads = cpu_count() - 2
    if num_threads < 1:
        num_threads = 1

    parser = argparse.ArgumentParser(
        description='Removing redundancy from the gb file.')
    parser.add_argument('gb', type=argparse.FileType('r'),
                   help='genbank file .gb')
    parser.add_argument('-n', '--num_threads', nargs='?', type = int,
                        const=num_threads, default=num_threads,
                    help="Number of threads to be executed in parallel.")
    
    return parser.parse_args()

def gb2fasta(gb):
    """Convert genbank file to fasta.

    This function converts a file in a genbank
    format to a file in a fasta file format.
    
    Parameters
    -------------
    gb: io.TextIOWrapper
        Genbank file

    Returns
    -------------
    fasta: str
        Name of the new fasta file

    """

    base = os.path.basename(gb.name)
    name = os.path.splitext(base)[0]

    fasta = "{}.fasta".format(name)
    SeqIO.convert(gb, "genbank", fasta, "fasta")

    return fasta

def clustering(fasta, num_threads):
    """Clustering sequences with VSEARCH

    This function uses VSEARCH suite to
    cluster the sequences using the cluster_fast
    algorithm.


    Parameters
    -------------
    fasta: str
        Fasta file name
    num_threads: int
        Number of threads to be used
    """

    vsearch = "vsearch -threads {num_threads} --cluster_fast {fasta}" \
            " --strand both --uc cluster.uc --id 1 --query_cov 1"
    subprocess.call(vsearch.format(num_threads = num_threads,
                                    fasta = fasta) , shell=True)

def removing_redundancy(gb):
    """Removing redundancy from genbank file after
    clustering sequences

    This function parses the uc file generated
    after clustering the sequences to discard
    the sequences from the original file.


    Parameters
    -------------
    gb: io.TextIOWrapper
        Genbank file
    """

    with open("cluster.uc") as cluster:    
        gb_all_data = SeqIO.index(gb.name, "genbank")

        base = os.path.basename(gb.name)
        name = os.path.splitext(base)[0]

        gb_c_name = "{}_c.gb".format(name)
        with open(gb_c_name, "w") as gb_c:
            for uc in tqdm(cluster, desc="Removing redundancy"):
                if uc[0] == "C":
                    row = uc.split("\t")
                    print(gb_all_data.get_raw(row[-2]).decode(),
                            file=gb_c, end="")


if __name__ == "__main__":
    args = getArguments()
    fasta = gb2fasta(args.gb)
    clustering(fasta, args.num_threads)
    print()
    removing_redundancy(args.gb)