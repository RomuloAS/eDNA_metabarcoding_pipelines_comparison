#!/usr/bin/env python3
import os
import glob
import argparse
import datetime
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from multiprocessing.dummy import Pool
from subprocess import run, check_output
from multiprocessing import Process, cpu_count

"""Remove adapters from fasta(q) files

From a list of files inside a folder,
it removes adapters from the 5', 3'
or both sides.
"""

DATETIME = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
PATH_RESULT = "{DATETIME}_Remove_Adapter".format(DATETIME=DATETIME)
CUTADAPT = "cutadapt {SIDE_ADAPTER}-o {OUTPUT} {INPUT}"

CHOICES_S = ["a", "b", "g"]

COMMANDS = {}

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
        description="Remove adapters from fasta(q) files")
    parser.add_argument("folder_fastaq", type=str,
                   help="A folder with fasta(q) files")
    parser.add_argument("folder_demultiplex", type=str,
                   help="A folder with tables with" \
                   " information about demultiplex adapters. "\
                   "Table files must have an identical name as fasta(q)")
    parser.add_argument("-n", "--num_threads", nargs="?", type = int,
                        const=num_threads, default=num_threads,
                    help="Number of threads to be executed in parallel.")
    parser.add_argument("-s", "--side", nargs="?", const="a",
                    default="a", choices=CHOICES_S,
                    type = lambda s : s.lower(),
                    help="Side of the sequence to cut the adapter"\
                    ", default: a (3' side) ")
    
    return parser.parse_args()

def parse_demultiplex_files(demultiplex_files):
    """Parse demultiplex files to get adapters.

    For each demultiplex table, parse the id of 
    the file and the adapter in normal and reverse
    direction
    
    demultiplex_files: list
        List with tables demultiplex info

    Returns
    -------------
    demultiplex: dataframe
        pandas dataframe with demultiplex adapters info
    """

    demultiplex = pd.DataFrame()

    for demultiplex_file in demultiplex_files:
        demultiplex = demultiplex.append(pd.read_csv("{DEMULTIPLEX_FILE}".format(
                                        DEMULTIPLEX_FILE = demultiplex_file), 
                                        sep = None, header = None,
                                        engine = 'python')
                                        )

    return demultiplex

def remove_adapters(fastaq_files, demultiplex, side, num_threads):
    """Remove adapters out of the fasta(q) files.

    From a list of fasta(q) files, this function
    removes adapters of the 5' end, 3' end
    or both sides.

    fastaq_files: list
        List of fasta(q) files
    demultiplex: dataframe
        Pandas dataframe with demultiplex info
    side: str
        Side of the sequence to remove adapters
    num_threads: int
        Number of threads for the multithreading
    """

    for fastaq_file in fastaq_files:
        file_name = Path(fastaq_file).stem
        base_name = file_name.split(".")[0]

        adapters = set(demultiplex.loc[demultiplex.iloc[:,1] == base_name,
                             2].values.flatten().tolist())

        side_adapter = ""
        for adapter in adapters:
            adapter = adapter.split(":")[0]

            if side == "a":
                side_adapter += "-{SIDE} {ADAPTER} ".format(
                    SIDE = side,
                    ADAPTER = str(Seq(adapter).reverse_complement()))

                continue

            side_adapter += "-{SIDE} {ADAPTER} ".format(
                    SIDE = side,
                    ADAPTER = adapter)

        COMMANDS["{BASE_NAME} R1".format(
                BASE_NAME = base_name)] = CUTADAPT.format(
                SIDE_ADAPTER = side_adapter,
                OUTPUT = "{PATH_RESULT}/Sequences/{BASE_NAME}"\
                    "_L001_R1_001.fastq.gz".format(
                    PATH_RESULT = PATH_RESULT,
                    BASE_NAME = base_name),
                INPUT = fastaq_file)

        COMMANDS["{BASE_NAME} R2".format(
                BASE_NAME = base_name)] = CUTADAPT.format(
                SIDE_ADAPTER = side_adapter,
                OUTPUT = "{PATH_RESULT}/Sequences/{BASE_NAME}"\
                    "_L001_R2_001.fastq.gz".format(
                    PATH_RESULT = PATH_RESULT,
                    BASE_NAME = base_name),
                INPUT = fastaq_file.replace("R1", "R2"))

    pool = Pool(num_threads)
    for returncode in pool.imap(execute_command, COMMANDS):
        if returncode:
            print("command failed: {}".format(returncode))


def execute_command(base_name):
    """Execute command

    This function executes command on terminal.

    base_name: str
        base_name of the fastq file
    """

    print(COMMANDS[base_name])
    with open("{PATH_RESULT}/Logs/{BASE_NAME}.log".format(
                PATH_RESULT = PATH_RESULT,
                BASE_NAME = base_name), "a") as log_file:
        run(COMMANDS[base_name], shell=True, stderr=log_file, stdout=log_file)

if __name__ == "__main__":
    args = getArguments()

    # Reading fastq files
    fastaq_files = sorted(glob.glob("{folder}/*R1*.fast*".format(
        folder=args.folder_fastaq)))
    demultiplex_files = sorted(glob.glob("{folder}/*.tsv".format(
        folder=args.folder_demultiplex)))
    side = args.side
    num_threads = args.num_threads

    # Parsing demultiplex files info
    demultiplex = parse_demultiplex_files(demultiplex_files)

    # Creating folders
    Path(PATH_RESULT).mkdir(parents=True, exist_ok=True)
    Path("{PATH_RESULT}/Logs".format(PATH_RESULT = PATH_RESULT)).mkdir(
            parents=True, exist_ok=True)
    Path("{PATH_RESULT}/Sequences".format(PATH_RESULT = PATH_RESULT)).mkdir(
            parents=True, exist_ok=True)

    # Removing adapters
    remove_adapters(fastaq_files, demultiplex, side, num_threads)