#!/usr/bin/env python3
import os
import sys
import pathlib
import argparse
import datetime
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from subprocess import check_output
from collections import defaultdict

"""Count reads from fasta or fastq inside the folder and
create a csv file.

For each fasta or fastq inside the folder informed count
the number of reads and create a table where the row name
is informed and column names are parsed from the files.
"""

DATETIME = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")

PATH_RESULT = "{DATETIME}_Count_Reads".format(DATETIME=DATETIME)

FASTA_COUNT_READS = '{Z}grep -c ">" {FASTA}'
FASTQ_COUNT_READS = "echo $({Z}cat {FASTQ}|wc -l)/4|bc"

CHOICES_F = ['fastq', 'fasta']

def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    parser = argparse.ArgumentParser(
        description="Count reads from FASTA or FASTQ files")
    parser.add_argument("folder_fastaq", type=str,
                   help="A folder with FASTA or FASTQ files")
    parser.add_argument("extension", type = lambda s : s.lower(), 
                    choices=CHOICES_F, help="Files extension")
    parser.add_argument("row_name", type=str,
                   help="Row name of the new table created")
    parser.add_argument("pattern", type=str,
                   help="A pattern to identify files to be parsed")
    parser.add_argument("-seqme", action='store_true',
                    help="Are the files from SeqME pipeline?")
    
    return parser.parse_args()

def count_reads(fastaq_files, row_name, pattern, extension, seqme):
    """Count the reads from FASTA/Q files and create a table

    From the list of FASTA or FASTQ files, this function
    counts the number of reads and create a csv table with 
    the row name informed and column names are parsed 
    from files.

    fastaq_files: list
        List with FASTA or FASTQ files
    row_name: str
        Row name
    pattern: str
        Pattern for the identification of files
    extension: str
        Extension of files
    seqme: bool
        Boolean if files are from SeqME
    """

    # Initialize dataframe
    df = pd.DataFrame()
    seqme_count = defaultdict(int)

    for fastaq_file in tqdm(fastaq_files, desc="Counting reads"):
        file_name = os.path.basename(fastaq_file)
        base_name = file_name.replace(pattern, "")

        if seqme:
            with open(fastaq_file) as fastaq:
                fasta = list(filter(None, fastaq.read().split(">")))

                for sequence in fasta:
                    base_name = sequence.split("|")[0]
                    seqme_count[base_name] = seqme_count[base_name] + 1

                continue

        z = ""
        if pathlib.Path(fastaq_file).suffix == ".gz":
            z = "z"

        try:
            if extension == "fasta":
                output = check_output(FASTA_COUNT_READS.format(
                                    Z=z,
                                    FASTA=fastaq_file),
                                    shell=True)
            elif extension == "fastq":
                output = check_output(FASTQ_COUNT_READS.format(
                                    Z=z,
                                    FASTQ=fastaq_file),
                                    shell=True)
        except:
            output = 0

        df.loc[row_name, base_name] = int(output)
    
    if seqme:
        df = pd.DataFrame(seqme_count, index=[row_name,])

    df.sort_index(axis=1).to_csv(
            "{PATH_RESULT}/{DATETIME}_Count_Reads.tsv".format(
            PATH_RESULT=PATH_RESULT,
            DATETIME=DATETIME),
                sep="\t")

if __name__ == "__main__":
    args = getArguments()

    Path(PATH_RESULT).mkdir(parents=True, exist_ok=True)

    fastaq_files = Path(args.folder_fastaq).rglob(
                        "*{PATTERN}*".format(PATTERN=args.pattern))

    count_reads(fastaq_files, args.row_name, args.pattern, args.extension, args.seqme)