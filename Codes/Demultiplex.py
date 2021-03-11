#!/usr/bin/env python3
import os
import sys
import csv
import glob
import shutil
import pathlib
import argparse
import datetime
import subprocess
from pathlib import Path

"""Demultiplex fastq files in a folder.

For each fastq file inside the folder,
collect the information from the table and
demultiplex it to new files.
"""

PATH_RESULT = "{}_Demultiplexed".format(datetime.datetime.now().strftime("%d%m%Y_%H%M%S"))

LOG = "{}/Demultiplex_error.log".format(PATH_RESULT)

DEMULTIPLEXING = "python demultiplex_obi_Sep_2017.py {TSV} {R1} {R2} {folder_to_save}"

def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    parser = argparse.ArgumentParser(
        description="Demultiplex FASTQ files"\
                    " using the information from the table")
    parser.add_argument("folder_fastq", type=str,
                   help="A folder with FASTQ files")
    parser.add_argument("folder_tables", type=str,
                   help="A folder with CSV/TSV table files")
    
    return parser.parse_args()

def demultiplexing(fastq_file, folder_tables):
    """Demultiplex FASTQ file.

    For each pair R1 R2, this function uses the table
    provided to demultiplex the fastq files.

    fastq_file: str
        fastq file name

    folder_tables: str
        Path to folder with csv/tsv files
    """

    file_name = os.path.basename(fastq_file)
    base_name = file_name.split("_")[0]

    table_name = "{dirname}/{name}*".format(dirname=folder_tables,
                                                    name=base_name)

    table = glob.glob(table_name)

    if not table:
        raise Exception("{table} table not found".format(
                            table=table))

    table = table[0]

    path_tsv = "{dirname}/TSV/".format(dirname=PATH_RESULT)
    Path(path_tsv).mkdir(parents=True, exist_ok=True)
    tsv_file = "{path_tsv}{table_name}.tsv".format(
                    path_tsv=path_tsv, table_name=base_name)

    with open(tsv_file, "w+") as output_tsv:
        with open(table) as input_table:
            table = input_table.read()
            sniffer = csv.Sniffer()
            dialect = sniffer.sniff(table)

            table = table.rstrip().lstrip().split("\n")
            csv.writer(output_tsv, delimiter='\t').writerows(
                csv.reader(table, delimiter=dialect.delimiter))


    path_demultiplexed = "{dirname}/{name}".format(
                            dirname=PATH_RESULT, name=base_name)

    Path(path_demultiplexed).mkdir(parents=True, exist_ok=True)

    subprocess.run(DEMULTIPLEXING.format(TSV=tsv_file, 
                                        R1=fastq_file,
                                        R2=fastq_file.replace("R1", "R2"),
                                        folder_to_save=path_demultiplexed),
                                        shell=True,
                                        stderr=sys.stderr,
                                        stdout=sys.stdout)


if __name__ == "__main__":
    args = getArguments()

    Path(PATH_RESULT).mkdir(parents=True, exist_ok=True)
    
    fastq_files = []
    fastq_files.extend(glob.glob("{folder}*R1*.fastq.gz".format(
        folder=args.folder_fastq)))

    for fastq_file in fastq_files:
        try:
            demultiplexing(fastq_file, args.folder_tables)
        except Exception as e:
            with open(LOG, "a") as log:
                print("{}\n".format(e), file=log)

    raw_reads = "{dirname}/raw_reads".format(dirname=PATH_RESULT)
    Path(raw_reads).mkdir(parents=True, exist_ok=True)

    print("\nMOVING to raw reads folder")
    fastq_demultiplexed = []
    fastq_demultiplexed.extend(sorted(glob.glob("{dirname}/**/*.fastq".format(
        dirname=PATH_RESULT), recursive=True)))

    fastq_files = {}
    for fastq_file in fastq_demultiplexed:
        file_name = os.path.basename(fastq_file)
        if file_name not in fastq_files:
            fastq_files[file_name] = [fastq_file]
        else:
            fastq_files[file_name].append(fastq_file)

    mv = "mv {fastq_file} {raw_reads}"
    cat = "cat {files} > {raw_reads}/{file_name}"
    for file_name, fastq_file in fastq_files.items():
        if len(fastq_files[file_name]) == 1:
            subprocess.run(mv.format(fastq_file=fastq_files[file_name][0], 
                                    raw_reads=raw_reads),
                                    shell=True, stderr=sys.stderr,
                                    stdout=sys.stdout)
        else:
            subprocess.run(cat.format(files=" ".join(fastq_files[file_name]), 
                                    raw_reads=raw_reads,
                                    file_name=file_name),
                                    shell=True, stderr=sys.stderr,
                                    stdout=sys.stdout)


    print("REMOVING invalid files and temporary folders")
    rm = "rm {raw_reads}/invalid*"
    subprocess.run(rm.format(raw_reads=raw_reads), shell=True,
                            stderr=sys.stderr, stdout=sys.stdout)

    folders = glob.glob("{PATH_RESULT}/*".format(PATH_RESULT=PATH_RESULT))
    for folder in folders:
        path = pathlib.PurePath(folder)
        if path.name != "raw_reads" \
                and path.name != "TSV":
            shutil.rmtree(folder)

    print("ZIPPING files")
    gzip = "gzip {raw_reads}/*.fastq"
    subprocess.run(gzip.format(raw_reads=raw_reads), shell=True,
                            stderr=sys.stderr, stdout=sys.stdout)