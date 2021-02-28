#!/usr/bin/env python3
import os
import argparse
import tempfile
import subprocess
from Bio import SeqIO
from tqdm import tqdm

"""
Trim alignment to remove large gaps in both extremities
of the sequences for building the phylogenetic tree.
The new trimmed alignment should not be used as a
reference to map reads.
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
        description="Trim alignment to remove large gaps in the" + 
                    " extremities of the sequences for building the tree.")
    parser.add_argument('aln', type=argparse.FileType('r'),
                   help='alignment file')
    
    return parser.parse_args()

def trim_alignment(aln):
    """Trim the alignment to remove large gaps

    This function uses trimal to remove large gaps from the
    multiple alignment sequences.

    Parameters
    -------------
    aln: io.TextIOWrapper
        Alignment file
    """

    base = os.path.basename(aln.name)
    name = os.path.splitext(base)[0]

    trimal = "trimal -in {} -gappyout > {}_t.aln".format(aln.name, name)
    subprocess.call(trimal , shell=True)

def convert2Phy(aln):
    """Convert Alignment to PHYLIP's format

    This function uses biopython package to read an alignment file,
    remove the name of the species (only accession number is kept),
    and convert it to PHYLIP's format.
    
    Parameters
    -------------
    aln: io.TextIOWrapper
        Alignment file
    """

    base = os.path.basename(aln.name)
    name = os.path.splitext(base)[0]

    aln_t = "{}_t.aln".format(name)
    phy = "{}.phy".format(name)

    with open(aln_t) as aln_file:
        aln_f = aln_file.read().split(">")

    for i in range(len(aln_f)):
        if aln_f[i]:
            aln_s = aln_f[i].split("_")
            aln_f[i] = aln_s[-1]
            if not aln_f[i][0].isalpha():
                aln_f[i] = "_".join(aln_s[-2:])

    with tempfile.NamedTemporaryFile(mode='w') as temp:
         temp.write(">".join(aln_f))
         temp.seek(0)

         SeqIO.convert(temp.name, "fasta", phy, "phylip-relaxed")



if __name__ == "__main__":
    args = getArguments()
    trim_alignment(args.aln)
    convert2Phy(args.aln)