#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import datetime
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from functools import partial
from multiprocessing.dummy import Pool
from multiprocessing import Process, cpu_count
from subprocess import run, check_output, DEVNULL, STDOUT

"""SeqME pipeline for metabarcoding detection.

For each sample fastq pair of files:
    
    .Joining paired-ends (fastq-join)
    .Quality filtering (fastx_toolkit)
    .Removing too short and too long sequences (Biopieces)
    .Clustering (USEARCH v10.0.240)
    .Creating an OTU table (USEARCH v10.0.240)
    .Alpha diversity & normalization (USEARCH v10.0.240)
    .Identifying OTUs by classifier (RDPTools)

"""

PATH_CLASSIFIER = "Classifier/rRNAClassifier.properties"

DATETIME = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")

PATH_RESULT = "{DATETIME}_SeqME".format(DATETIME=DATETIME)
PATH_FASTQ_JOINED = "{PATH_RESULT}/Fastq_Joined"
PATH_QUALITY_FILTERED = "{PATH_RESULT}/Fastq_Quality_Filtered"
PATH_FASTA = "{PATH_RESULT}/Fasta"
PATH_FASTA_REMOVED_SHORT_LONG_SEQ = "{PATH_RESULT}/Fasta_Removed_Short_Long_Seq"
PATH_UNIQUES = "{PATH_RESULT}/Fasta_Uniques"
PATH_OTUS_FASTA = "{PATH_RESULT}/Fasta_OTUs"
PATH_OTUS_TABLE = "{PATH_RESULT}/Table_OTUs"
PATH_OTUS_TABLE_NORMALIZE = "{PATH_RESULT}/Table_OTUs_Normalized"
PATH_ALPHA_DIVERSITY = "{PATH_RESULT}/Table_Alpha_Diversity"
PATH_OTU_IDENTIFIED = "{PATH_RESULT}/Table_OTUs_identified"

FASTQ_JOIN = "fastq-join -v ' ' -p 15 -m 15 {R1} {R2} -o {FASTQ_OUTPUT}"
FASTQ_QUALITY_FILTER = "fastq_quality_filter -i {FASTQ_INPUT}"\
                        " -Q33 -q 20 -p 50 -o {FASTQ_OUTPUT}"
FASTQ_TO_FASTA = "fastq_to_fasta -i {FASTQ_INPUT} -o {FASTA_OUTPUT}"
FASTA_REMOVE_SHORT_LONG_SEQ = "read_fasta -i {FASTA_INPUT}" \
                                " | grab -e 'SEQ_LEN >= 90'"\
                                " | grab -e 'SEQ_LEN <= 150'"\
                                " | write_fasta -x -o {FASTA_OUTPUT}"
FASTA_CLUSTER_UNIQUES = "usearch -fastx_uniques {FASTA_INPUT}"\
            " -fastaout {FASTA_OUTPUT} -uc {UC_OUTPUT}"\
            " -sizeout -relabel Uniq"
FASTA_CLUSTER_OTU = "usearch -cluster_otus {FASTA_INPUT}"\
                    " -otus {FASTA_OUTPUT} -relabel Otu"
FASTA_CREATE_OTU_TABLE = "usearch -otutab {FASTA_INPUT}"\
                        " -otus {FASTA_OTU_INPUT} -otutabout"\
                        " {OTU_OUTPUT} -mapout {MAP_OUTPUT}"
OTU_NORMALIZE = "usearch -otutab_rare {OTU_INPUT} -sample_size"\
                " 5000 -output {OTU_OUTPUT}"
ALPHA_DIVERSITY = "usearch -alpha_div {OTU_INPUT}"\
                    " -output {ALPHA_OUTPUT}"
OTU_IDENTIFY = "classifier classify -t {CLASSIFIER}"\
                " -c 1 -w 150 -o {TABLE_OUTPUT}"\
                " -h {HIER_OUTPUT} {FASTA_INPUT}"

GUNZIP = "gunzip {folder}/*.gz"

COUNT_READS = "echo $(cat {fastq}|wc -l)/4|bc"

THRESHOLDS = np.round(np.arange(0.1, 1.01, 0.01), 2)

CHOICES_R = ["superkingdom", "kingdom", "phylum",
             "class", "order", "family", "genus",
             "species"]


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
        description="SeqME pipeline for Metabarcoding")
    parser.add_argument("folder_fastq", type=str,
                   help="A folder with fastq files")
    parser.add_argument("-t", "--threshold", nargs="?", const=1, default=1,
                    choices=THRESHOLDS, type = float,
                    help="Specify the minimum threshold"\
                    " to the taxonomy rank be kept, default: 1")
    parser.add_argument("-o", "--only_joining", nargs="?", const="",
                        default="", type = str,
                    help="Inform folder - Only the final joining of"\
                    " the results is done")
    parser.add_argument("-n", "--num_threads", nargs="?", type = int,
                        const=num_threads, default=num_threads,
                    help="Number of threads to be executed in parallel")
    parser.add_argument("-no", "--normalized", action="store_true",
                    help="Use normalized data")
    parser.add_argument("-r", "--rank", nargs="?", const="species", default="species",
                    choices=CHOICES_R, type = lambda s : s.lower(),
                    help="Lowest taxonomic classification rank"\
                    " to be in the result, default: species")
    
    return parser.parse_args()

def join_paired_ends(fastq_file, base_name):
    """Join paired-end Illumina data.

    From a pair forward (R1) and reverse (R2),
    this function creates the command line 
    to merge the pair of files into a single
    sequence using fastq-join.
    
    fastq_file: str
        fastq file
    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output file name
    output_file_name = "{}_\%.fastq".format(base_name)

    # Folder to be saved
    fastq_output = "{PATH_FASTQ_JOINED}/{FASTQ_OUTPUT}".format(
                    PATH_FASTQ_JOINED=PATH_FASTQ_JOINED,
                    FASTQ_OUTPUT= output_file_name
                    )

    return ("Joining paired ends: {file_name}".format(file_name=base_name),
            "{PATH_FASTQ_JOINED}/{base_name}.log".format(
                PATH_FASTQ_JOINED=PATH_FASTQ_JOINED,
                base_name=base_name),
            FASTQ_JOIN.format(R1=fastq_file, 
                                R2=fastq_file.replace("R1", "R2"),
                                FASTQ_OUTPUT=fastq_output)
        )

def quality_filtering(base_name):
    """Quality filtering of fastq file.

    From a fastq file, this function creates the 
    command line to remove low quality nucleotides.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and joined file names
    output_file_name = "{}.fastq".format(base_name)
    joined_file_name = "{}_join.fastq".format(base_name)

    # Joined file path
    joined_file_path = "{PATH_FASTQ_JOINED}/{FASTQ_INPUT}".format(
                    PATH_FASTQ_JOINED=PATH_FASTQ_JOINED,
                    FASTQ_INPUT= joined_file_name
                    )

    # Folder to be saved
    fastq_output = "{PATH_QUALITY_FILTERED}/{FASTQ_OUTPUT}".format(
                    PATH_QUALITY_FILTERED=PATH_QUALITY_FILTERED,
                    FASTQ_OUTPUT= output_file_name
                    )

    return ("Quality filtering: {file_name}".format(file_name=base_name),
            "{PATH_QUALITY_FILTERED}/{base_name}.log".format(
                PATH_QUALITY_FILTERED=PATH_QUALITY_FILTERED,
                base_name=base_name),
            FASTQ_QUALITY_FILTER.format(FASTQ_INPUT=joined_file_path,
                                FASTQ_OUTPUT=fastq_output)
        )

def convert_fastq_to_fasta(base_name):
    """Convert fastq to fasta.

    From a fastq file, this function converts fastq
    to fasta format.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    input_file_name = "{}.fastq".format(base_name)
    output_file_name = "{}.fasta".format(base_name)

    fastq_input = "{PATH_QUALITY_FILTERED}/{FASTQ_INPUT}".format(
                    PATH_QUALITY_FILTERED=PATH_QUALITY_FILTERED,
                    FASTQ_INPUT= input_file_name
                    )
    fasta_output = "{PATH_FASTA}/{FASTA_OUTPUT}".format(
                    PATH_FASTA=PATH_FASTA,
                    FASTA_OUTPUT= output_file_name
                    )

    return ("Converting fastq to fasta: {file_name}".format(file_name=base_name),
            "{PATH_FASTA}/{base_name}.log".format(
                PATH_FASTA=PATH_FASTA,
                base_name=base_name),
            FASTQ_TO_FASTA.format(FASTQ_INPUT=fastq_input,
                                FASTA_OUTPUT=fasta_output)
        )

def remove_short_long_seq(base_name):
    """Remove too short and too long sequences.

    From a fasta file, this function removes
    too short and too long sequences.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    fasta_file_name = "{}.fasta".format(base_name)

    fasta_input = "{PATH_FASTA}/{FASTA_INPUT}".format(
                PATH_FASTA=PATH_FASTA,
                FASTA_INPUT= fasta_file_name
                    )
    fasta_output = "{PATH_FASTA_REMOVED_SHORT_LONG_SEQ}/{FASTA_OUTPUT}".format(
                PATH_FASTA_REMOVED_SHORT_LONG_SEQ=PATH_FASTA_REMOVED_SHORT_LONG_SEQ,
                FASTA_OUTPUT= fasta_file_name
                    )

    return ("Removing short and long seq: {file_name}".format(file_name=base_name),
            "{PATH_FASTA_REMOVED_SHORT_LONG_SEQ}/{base_name}.log".format(
                PATH_FASTA_REMOVED_SHORT_LONG_SEQ=PATH_FASTA_REMOVED_SHORT_LONG_SEQ,
                base_name=base_name),
            FASTA_REMOVE_SHORT_LONG_SEQ.format(FASTA_INPUT=fasta_input,
                                                FASTA_OUTPUT=fasta_output)
        )

def cluster_uniques(base_name):
    """Cluster uniques sequences (dereplication).

    From a fasta file, this function finds a set
    of unique sequences in the file.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    fasta_file_name = "{}.fasta".format(base_name)
    uc_file_name = "{}.uc".format(base_name)

    fasta_input = "{PATH_FASTA_REMOVED_SHORT_LONG_SEQ}/{FASTA_INPUT}".format(
                PATH_FASTA_REMOVED_SHORT_LONG_SEQ=PATH_FASTA_REMOVED_SHORT_LONG_SEQ,
                FASTA_INPUT= fasta_file_name
                    )
    fasta_output = "{PATH_UNIQUES}/{FASTA_OUTPUT}".format(
                PATH_UNIQUES=PATH_UNIQUES,
                FASTA_OUTPUT= fasta_file_name
                    )
    uc_output = "{PATH_UNIQUES}/{UC_OUTPUT}".format(
                PATH_UNIQUES=PATH_UNIQUES,
                UC_OUTPUT= uc_file_name
                    )

    return ("Clustering uniques (Dereplication): {file_name}".format(file_name=base_name),
            "{PATH_UNIQUES}/{base_name}.log".format(
                PATH_UNIQUES=PATH_UNIQUES,
                base_name=base_name),
            FASTA_CLUSTER_UNIQUES.format(FASTA_INPUT=fasta_input,
                                        FASTA_OUTPUT=fasta_output,
                                        UC_OUTPUT=uc_output)
        )

def cluster_otus(base_name):
    """Cluster OTUs.

    From a fasta file, this function does a OTU
    clustering. Chimeras are also filtered
    during this step.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    fasta_file_name = "{}.fasta".format(base_name)

    fasta_input = "{PATH_UNIQUES}/{FASTA_INPUT}".format(
                PATH_UNIQUES=PATH_UNIQUES,
                FASTA_INPUT= fasta_file_name
                    )
    fasta_output = "{PATH_OTUS_FASTA}/{FASTA_OUTPUT}".format(
                PATH_OTUS_FASTA=PATH_OTUS_FASTA,
                FASTA_OUTPUT= fasta_file_name
                    )

    return ("Clustering OTUs: {file_name}".format(file_name=base_name),
            "{PATH_OTUS_FASTA}/{base_name}.log".format(
                PATH_OTUS_FASTA=PATH_OTUS_FASTA,
                base_name=base_name),
            FASTA_CLUSTER_OTU.format(FASTA_INPUT=fasta_input,
                                    FASTA_OUTPUT=fasta_output)
        )

def create_otu_tables(base_name):
    """Create OTU tables.

    From a fasta file, this function creates
    an OTU table with the identification of which
    OTU the sequence belongs to, and number of sequences
    for each OTU.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    fasta_file_name = "{}.fasta".format(base_name)
    otu_file_name = "{}.otu".format(base_name)
    map_file_name = "{}.map".format(base_name)

    fasta_input = "{PATH_FASTA_REMOVED_SHORT_LONG_SEQ}/{FASTA_INPUT}".format(
                PATH_FASTA_REMOVED_SHORT_LONG_SEQ=PATH_FASTA_REMOVED_SHORT_LONG_SEQ,
                FASTA_INPUT= fasta_file_name
                    )
    fasta_otu_input = "{PATH_OTUS_FASTA}/{FASTA_OTU_INPUT}".format(
                PATH_OTUS_FASTA=PATH_OTUS_FASTA,
                FASTA_OTU_INPUT= fasta_file_name
                    )
    otu_output = "{PATH_OTUS_TABLE}/{OTU_OUTPUT}".format(
                PATH_OTUS_TABLE=PATH_OTUS_TABLE,
                OTU_OUTPUT= otu_file_name
                    )
    map_output = "{PATH_OTUS_TABLE}/{MAP_OUTPUT}".format(
                PATH_OTUS_TABLE=PATH_OTUS_TABLE,
                MAP_OUTPUT= map_file_name
                    )

    return ("Creating OTU tables: {file_name}".format(file_name=base_name),
            "{PATH_OTUS_TABLE}/{base_name}.log".format(
                PATH_OTUS_TABLE=PATH_OTUS_TABLE,
                base_name=base_name),
            FASTA_CREATE_OTU_TABLE.format(FASTA_INPUT=fasta_input,
                                        FASTA_OTU_INPUT=fasta_otu_input,
                                        OTU_OUTPUT=otu_output,
                                        MAP_OUTPUT=map_output)
        )

def normalize_otu_tables(base_name):
    """Normalize OTU tables.

    From a fasta file, this function
    normalizes all samples.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    otu_file_name = "{}.otu".format(base_name)

    otu_input = "{PATH_OTUS_TABLE}/{OTU_INPUT}".format(
                PATH_OTUS_TABLE=PATH_OTUS_TABLE,
                OTU_INPUT= otu_file_name
                    )
    otu_output = "{PATH_OTUS_TABLE_NORMALIZE}/{OTU_OUTPUT}".format(
                PATH_OTUS_TABLE_NORMALIZE=PATH_OTUS_TABLE_NORMALIZE,
                OTU_OUTPUT= otu_file_name
                    )

    return ("Normalizing OTU tables: {file_name}".format(file_name=base_name),
            "{PATH_OTUS_TABLE_NORMALIZE}/{base_name}.log".format(
                PATH_OTUS_TABLE_NORMALIZE=PATH_OTUS_TABLE_NORMALIZE,
                base_name=base_name),
            OTU_NORMALIZE.format(OTU_INPUT=otu_input,
                                OTU_OUTPUT=otu_output)
        )

def calculate_alpha_diversity(base_name):
    """Calculate alpha diversity.

    From a fasta file, this function calculates
    the alpha diversity.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    otu_file_name = "{}.otu".format(base_name)
    alpha_file_name = "{}.alpha".format(base_name)

    otu_input = "{PATH_OTUS_TABLE_NORMALIZE}/{OTU_INPUT}".format(
                PATH_OTUS_TABLE_NORMALIZE=PATH_OTUS_TABLE_NORMALIZE,
                OTU_INPUT= otu_file_name
                    )
    alpha_output = "{PATH_ALPHA_DIVERSITY}/{ALPHA_OUTPUT}".format(
                PATH_ALPHA_DIVERSITY=PATH_ALPHA_DIVERSITY,
                ALPHA_OUTPUT= alpha_file_name
                    )

    return ("Calculating alpha diversity: {file_name}".format(file_name=base_name),
            "{PATH_ALPHA_DIVERSITY}/{base_name}.log".format(
                PATH_ALPHA_DIVERSITY=PATH_ALPHA_DIVERSITY,
                base_name=base_name),
            ALPHA_DIVERSITY.format(OTU_INPUT=otu_input,
                                ALPHA_OUTPUT=alpha_output)
        )

def identify_otu_by_classifier(base_name):
    """Identify OTU by classifier.

    From a fasta file, this function identifies
    taxonomy rank for each OTU using a Naive Bayes
    classifier.

    base_name: str
        fastq file base name

    Returns
    -------------
    command: tuple
        A tuple with name, command and log
    """

    # Output and input file names
    fasta_file_name = "{}.fasta".format(base_name)
    hier_file_name = "{}.hier".format(base_name)
    table_file_name = "{}.csv".format(base_name)

    fasta_otu_input = "{PATH_OTUS_FASTA}/{FASTA_OTU_INPUT}".format(
                PATH_OTUS_FASTA=PATH_OTUS_FASTA,
                FASTA_OTU_INPUT= fasta_file_name
                    )
    hier_output = "{PATH_OTU_IDENTIFIED}/{HIER_OUTPUT}".format(
                PATH_OTU_IDENTIFIED=PATH_OTU_IDENTIFIED,
                HIER_OUTPUT= hier_file_name
                    )
    table_output = "{PATH_OTU_IDENTIFIED}/{TABLE_OUTPUT}".format(
                PATH_OTU_IDENTIFIED=PATH_OTU_IDENTIFIED,
                TABLE_OUTPUT= table_file_name
                    )

    return ("Identifying OTU by classifier: {file_name}".format(file_name=base_name),
            "{PATH_OTU_IDENTIFIED}/{base_name}.log".format(
                PATH_OTU_IDENTIFIED=PATH_OTU_IDENTIFIED,
                base_name=base_name),
            OTU_IDENTIFY.format(CLASSIFIER=PATH_CLASSIFIER,
                                TABLE_OUTPUT=table_output,
                                HIER_OUTPUT=hier_output,
                                FASTA_INPUT=fasta_otu_input)
        )

def join_results(base_names, threshold, normalized, rank_to_stop):
    """Join the results in one unique csv file.

    From the list of fastq files, this function joins
    them and create a unique file with the taxonomic rank
    and the number of reads.

    base_names: list
        List with base names of fastq files
    threshold: float
        Minimum threshold to the taxonomy rank be kept
    normalized: bool
        Use normalized data
    rank_to_stop: str
        Taxonomic rank to stop searching
    """

    # Initialize dataframe
    df = pd.DataFrame()

    otus_identified = sorted(glob.glob(
        "{PATH_OTU_IDENTIFIED}/*.csv".format(
        PATH_OTU_IDENTIFIED=PATH_OTU_IDENTIFIED)))

    for otu_identified in tqdm(otus_identified, desc="Joining results in a unique file"):
        try:
            file_name = os.path.basename(otu_identified)
            base_name = os.path.splitext(file_name)[0]

            with open("{PATH_OTU_IDENTIFIED}/{TABLE}.csv".format(
                    PATH_OTU_IDENTIFIED=PATH_OTU_IDENTIFIED,
                    TABLE=base_name
                        )) as table_file:
                table = table_file.readlines()

            if not normalized:
                otu_path = "{PATH_OTUS_TABLE}/{OTU}.otu".format(
                        PATH_OTUS_TABLE=PATH_OTUS_TABLE,
                        OTU=base_name
                            )
            else:
                otu_path = "{PATH_OTUS_TABLE_NORMALIZE}/{OTU}.otu".format(
                        PATH_OTUS_TABLE_NORMALIZE=PATH_OTUS_TABLE_NORMALIZE,
                        OTU=base_name
                            )

            with open(otu_path) as otu_file:
                df = df.append(pd.read_table(otu_file, index_col=0,
                                        names=[base_name], skiprows=1))

            for line in table:
                line = line.split("\t")
                otu_id = line[0]
                
                for i in range(len(line)-1, 1, -3):
                    value = float(line[i])
                    rank = line[i-1]
                    tax = line[i-2]

                    if value >= threshold:
                        df.rename(index={otu_id:tax}, inplace=True)
                        break

                    if rank == rank_to_stop.lower():
                        break

                if otu_id in df.index:
                    df.drop(otu_id, inplace=True)

        except Exception as e:
            with open("{PATH_RESULT}/{DATETIME}_join_results.log".format(
                PATH_RESULT=PATH_RESULT,
                DATETIME=DATETIME),
                "a") as log:
                print(e, file=log)

    df.index.name = "TAX"
    df = df.groupby(df.index).sum()
    df.sort_index(axis=1, inplace=True)
    df.sort_index().to_csv("{PATH_RESULT}/{DATETIME}_SeqME.tsv".format(
            PATH_RESULT=PATH_RESULT,
            DATETIME=DATETIME),
                sep="\t")

def initialize_paths():
    """Initialize Paths

    This function initilizes paths to
    the results of each step of the pipeline.
    """

    global PATH_FASTQ_JOINED, PATH_QUALITY_FILTERED, PATH_FASTA, \
            PATH_FASTA_REMOVED_SHORT_LONG_SEQ, PATH_UNIQUES, \
            PATH_OTUS_FASTA, PATH_OTUS_TABLE, PATH_OTUS_TABLE_NORMALIZE, \
            PATH_ALPHA_DIVERSITY, PATH_OTU_IDENTIFIED

    PATH_FASTQ_JOINED = PATH_FASTQ_JOINED.format(PATH_RESULT=PATH_RESULT)
    PATH_QUALITY_FILTERED = PATH_QUALITY_FILTERED.format(PATH_RESULT=PATH_RESULT)
    PATH_FASTA = PATH_FASTA.format(PATH_RESULT=PATH_RESULT)
    PATH_FASTA_REMOVED_SHORT_LONG_SEQ = PATH_FASTA_REMOVED_SHORT_LONG_SEQ.format(
                                        PATH_RESULT=PATH_RESULT)
    PATH_UNIQUES = PATH_UNIQUES.format(PATH_RESULT=PATH_RESULT)
    PATH_OTUS_FASTA = PATH_OTUS_FASTA.format(PATH_RESULT=PATH_RESULT)
    PATH_OTUS_TABLE = PATH_OTUS_TABLE.format(PATH_RESULT=PATH_RESULT)
    PATH_OTUS_TABLE_NORMALIZE = PATH_OTUS_TABLE_NORMALIZE.format(
                                PATH_RESULT=PATH_RESULT)
    PATH_ALPHA_DIVERSITY = PATH_ALPHA_DIVERSITY.format(
                            PATH_RESULT=PATH_RESULT)
    PATH_OTU_IDENTIFIED = PATH_OTU_IDENTIFIED.format(PATH_RESULT=PATH_RESULT)

def create_folders():
    """Create folders

    This function create folders to
    the results of each step of the pipeline
    """

    Path(PATH_RESULT).mkdir(parents=True, exist_ok=True)
    Path(PATH_FASTQ_JOINED).mkdir(parents=True, exist_ok=True)
    Path(PATH_QUALITY_FILTERED).mkdir(parents=True, exist_ok=True)
    Path(PATH_FASTA).mkdir(parents=True, exist_ok=True)
    Path(PATH_FASTA_REMOVED_SHORT_LONG_SEQ).mkdir(parents=True, exist_ok=True)
    Path(PATH_UNIQUES).mkdir(parents=True, exist_ok=True)
    Path(PATH_OTUS_FASTA).mkdir(parents=True, exist_ok=True)
    Path(PATH_OTUS_TABLE).mkdir(parents=True, exist_ok=True)
    Path(PATH_OTUS_TABLE_NORMALIZE).mkdir(parents=True, exist_ok=True)
    Path(PATH_ALPHA_DIVERSITY).mkdir(parents=True, exist_ok=True)
    Path(PATH_OTU_IDENTIFIED).mkdir(parents=True, exist_ok=True)

def execute_command(base_name):
    """Execute command

    This function executes command from terminal.

    base_name: str
        Base name of the fastq file
    """

    for command in commands[base_name]:
        print(command[0])
        with open(command[1], "a") as log_file:
            run(command[2], shell=True, stderr=log_file, stdout=log_file)


if __name__ == "__main__":
    args = getArguments()

    initial_time = datetime.datetime.now()

    # Unzipping compressed files
    run(GUNZIP.format(folder=args.folder_fastq),
        shell=True, stdout=DEVNULL, stderr=STDOUT)

    # Reading fastq files
    fastq_files = sorted(glob.glob("{folder}/*R1*.fastq".format(
        folder=args.folder_fastq)))
    base_names = []
    commands = {}
    reads = {}

    if args.only_joining:
        PATH_RESULT = args.only_joining
        initialize_paths()
    else:
        initialize_paths()

    for fastq_file in tqdm(fastq_files, desc="Reading fastq files"):
        file_name = os.path.basename(fastq_file)
        base_name = file_name.split("_")[0].split(".R1")[0]
        base_names.append(base_name)

        output = check_output(COUNT_READS.format(
                                fastq=fastq_file), shell=True)
        reads[base_name] = int(output)

        commands[base_name] = [join_paired_ends(fastq_file, base_name),
                                quality_filtering(base_name),
                                convert_fastq_to_fasta(base_name),
                                remove_short_long_seq(base_name),
                                cluster_uniques(base_name),
                                cluster_otus(base_name),
                                create_otu_tables(base_name),
                                normalize_otu_tables(base_name),
                                calculate_alpha_diversity(base_name),
                                identify_otu_by_classifier(base_name)]

    if args.only_joining:
        if Path(args.only_joining).exists():
            join_results(base_names, args.threshold, args.normalized, args.rank)
        else:
            print("Folder '{FOLDER}' does not exist".format(
                FOLDER=args.only_joining))

        raise SystemExit(0)

    # Creating folders
    create_folders()

    # Running commands
    num_threads = args.num_threads

    pool = Pool(num_threads)
    for returncode in pool.imap(execute_command, commands):
        if returncode:
            print("command failed: {}".format(returncode))
    

    join_results(base_names, args.threshold, args.normalized, args.rank)

    final_time = datetime.datetime.now()
    total_time = final_time - initial_time
    print("\nTotal time to run the pipeline: {TOTAL_TIME}\n".format(
            TOTAL_TIME=total_time))
