#!/usr/bin/env python3
import os
import pathlib
import argparse
import tempfile
import datetime
from tqdm import tqdm
from Bio import SeqIO
from Bio import Entrez
from ete3 import NCBITaxa
from Bio.Align.Applications import *
from multiprocessing import cpu_count

"""
Sequences Alignment: 

Align sequences by group of sequences
based on taxonomic rank chosen (example: species).

A file with primers can be provided to be included
in the alignment.

Multiple sequence alignment program can be chosen.
"""

CHOICES_P = ["mafft", "clustalo", "muscle", "all"]

CHOICES_R = ["superkingdom", "kingdom", "phylum",
            "subphylum", "superclass", "class",
            "subclass", "order", "suborder",
            "family", "subfamily", "tribe",
            "genus", "species", "all"]

CHOICES_S = ["taxdump", "ncbi"]

LINEAGES = {}

PATH_MAIN = "{}_alignments/".format(datetime.datetime.now().strftime("%d%m%Y_%H%M%S"))

LOG = "{}Alignment.log".format(PATH_MAIN)


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
        description="Alignment of sequences from a genbank"\
                    " file based on taxonomic rank.")
    parser.add_argument("gb", type=argparse.FileType("r"),
                   help="genbank file (.gb)")
    parser.add_argument("-p", "--program", nargs="*", default=["mafft"],
                    choices=CHOICES_P, type = lambda s : s.lower(),
                    help="Multiple sequence alignment program, default: mafft")
    parser.add_argument("-r", "--rank", nargs="*", default=["species"],
                    choices=CHOICES_R, type = lambda s : s.lower(),
                    help="Taxonomic classification rank to be used"\
                    " to separate the groups, default: species")
    parser.add_argument("-s", "--source", nargs="?", const="taxdump", default="taxdump",
                    choices=CHOICES_S, type = lambda s : s.lower(),
                    help="Source to be used to collect"\
                    " the info about the taxonomic rank, default: taxdump")
    parser.add_argument("-pr", "--primers", type=argparse.FileType("r"),
                    help="A fasta file with primers")
    parser.add_argument("-sp", "--species_from_file", action="store_true",
                    help="Should the species name from file"\
                    " be used or be collected from NCBI/taxdump?")
    parser.add_argument('-n', '--num_threads', nargs='?', type = int,
                    const=num_threads, default=num_threads,
                    help="Number of threads to be executed in parallel.")
    
    return parser.parse_args()

def get_tax_lineage(taxonid, source):
    """Return taxonomy lineage information

    This function uses Biopython library to connect NCBI database
    and search for taxonomy information or ete3 to download
    taxdump file and search the information locally.

    Parameters
    -------------
    taxonid : string
        Taxonomic id of the species
    source : string
        Source to be used to collect the info about the taxonid

    Returns
    -------------
    lineage: dict
        Species lineage

    """

    if taxonid not in LINEAGES:
        if source == "taxdump":
            ncbi_taxdump = NCBITaxa()
            lineage_ids = ncbi_taxdump.get_lineage(taxonid)
            ranks = ncbi_taxdump.get_rank(lineage_ids)
            names = ncbi_taxdump.get_taxid_translator(lineage_ids)
            lineage = {ranks[i]:names[i] for i in lineage_ids}

            LINEAGES[taxonid] = lineage
            return LINEAGES[taxonid]

        while True:
            data = ""
            try:
                Entrez.email = "Your.Name.Here@example.org"
                handle = Entrez.efetch(id = taxonid, db = "taxonomy", retmode = "xml")    
                data = Entrez.read(handle)
                handle.close()
            except Exception as e:
                with open(LOG, "a") as log:
                    print("Error when searching information about {}".format(taxonid),
                        file=log)

            if data:
                break
    
        lineage = {d["Rank"]:d["ScientificName"] for d in data[0]["LineageEx"]}
        lineage[data[0]["Rank"]] = data[0]["ScientificName"]
        LINEAGES[taxonid] = lineage

    
    return LINEAGES[taxonid]

def read_sequences(gb, rank, source, species_from_file):
    """Read the genbank file and parse the sequences
    for the taxonomic rank

    This function uses Biopython library to scan the genbank file
    and parse the sequences for the taxonomic rank.

    Parameters
    -------------
    gb : io.TextIOWrapper
        A genbank file
    rank: string
        Taxonomic rank
    source : string
        Source to be used to collect the info about the taxonid
    species_from_file: bool
        Indicate if species name from file should be used

    Returns
    -------------
    sequences: dictionary
        A dictionary with key representing ranks and values
        representing sequences in fasta format with species name
        as the header of the sequence

    """

    print("Rank: {}".format(rank))
    sequences = {}
    for r in tqdm(gb, desc="Reading sequences"):
        record = gb.get(r)
        for feature in record.features:
            if feature.type == "source" and \
                "taxon" in feature.qualifiers["db_xref"][0]:
                taxonid = feature.qualifiers["db_xref"][0].split(":")[1]

        lineage = get_tax_lineage(taxonid, source)

        if species_from_file:
            lineage["species"] = record.features[0].qualifiers["organism"][0]
        else:
            if "species" not in lineage:
                lineage["species"] = record.features[0].qualifiers["organism"][0]

        try:
            if lineage[rank] not in sequences:
                sequences[lineage[rank]] = [">{}_{}\n{}".format(
                                                lineage["species"].replace(" ", "_"),
                                                record.id, record.seq)]
            else:
                sequences[lineage[rank]].append(">{}_{}\n{}".format(
                                                lineage["species"].replace(" ", "_"),
                                                record.id, record.seq))
        except:
            with open(LOG, "a") as log:
                print("\nRank '{}' not found for organism '{}', taxonid '{}'".format(
                    rank, lineage["species"], taxonid), file=log)
            
    return sequences

def alignment(sequences, program, rank, primers, num_threads):
    """Sequences alignment using the program chosen by rank level.

    This function uses either MAFFT, Clustal omega, or Muscle
    to perform a multiple alignment for each group
    of sequences in the rank(s) chosen.

    Parameters
    -------------
    sequences : dictionary
        A dictionary with sequences
    program: string
        Program to be used to align the sequences
    rank: string
        Taxonomic rank
    primers: list
        A list of primers to be included in the alignment
    num_threads: int
        Number of threads to be used

    """

    print("Program: {}".format(program))
    path_alignments = "{}{}_alignments/{}".format(PATH_MAIN, program, rank)
    pathlib.Path(path_alignments).mkdir(parents=True, exist_ok=True)

    for seq in tqdm(sequences, desc="Alignment"):
        with tempfile.NamedTemporaryFile(mode="w") as temp:
            temp.write("\n".join(primers + sorted(sequences[seq])))
            temp.seek(0)

            path_seq = "{}/{}".format(path_alignments, seq.replace(" ", "_"))
            pathlib.Path(path_seq).mkdir(parents=True, exist_ok=True)

            try:
                if program == "clustalo":
                    cmdline = ClustalOmegaCommandline(
                                        infile=temp.name,
                                        outfile="{}/{}.aln".format(path_seq,
                                        seq.replace(" ", "_")),
                                        guidetree_out="{}/{}.tree".format(path_seq,
                                        seq.replace(" ", "_")),
                                        force=True,
                                        threads = num_threads
                                        )
                    cmdline()

                elif program == "mafft":
                    cmdline = MafftCommandline(input=temp.name, treeout=True,
                                                localpair=True, maxiterate=1000,
                                                adjustdirectionaccurately=True,
                                                thread = num_threads)
                    stdout, stderr = cmdline()
                    os.rename("{}.tree".format(temp.name), "{}/{}.tree".format(path_seq,
                                            seq.replace(" ", "_")))
                    with open("{}/{}.aln".format(path_seq,
                                            seq.replace(" ", "_")), "w") as align:
                            print(stdout, file=align)

                elif program == "muscle":
                    cmdline = MuscleCommandline(input=temp.name,
                                                clwout="{}/{}.clw".format(path_seq,
                                                seq.replace(" ", "_")),
                                                fastaout="{}/{}.fa".format(path_seq,
                                                seq.replace(" ", "_")),
                                                htmlout="{}/{}.html".format(path_seq,
                                                seq.replace(" ", "_")),
                                                msfout="{}/{}.msf".format(path_seq,
                                                seq.replace(" ", "_")),
                                                phyiout="{}/{}.phy".format(path_seq,
                                                seq.replace(" ", "_")),
                                                tree1="{}/{}_1.tree".format(path_seq,
                                                seq.replace(" ", "_")),
                                                tree2="{}/{}_2.tree".format(path_seq,
                                                seq.replace(" ", "_"))
                                                )

                    cmdline()

            except Exception as e:
                if os.listdir(path_seq) == []:
                    os.rmdir(path_seq)
                with open(LOG, "a") as log:
                    print("{species}: {error}".format(
                            species=seq, error=e), file=log)
                continue


if __name__ == "__main__":
    args = getArguments()
    
    pathlib.Path(PATH_MAIN).mkdir(parents=True, exist_ok=True)
    gb = SeqIO.index(args.gb.name, "genbank")
    with open(LOG, "w"): pass

    primers = []
    if args.primers:
        try:
            primers = [">{}\n{}".format(r.description.replace(" ", "_"), r.seq)
                        for r in SeqIO.parse(args.primers, "fasta")]
        except:
            with open(LOG, "a") as log:
                print("Fasta primers not found!! " +
                        "Continuing without primers", file=log)

    ranks = args.rank
    if "all" in args.rank:
        ranks = CHOICES_R[:-1]

    programs = args.program
    if "all" in args.program:
        programs = CHOICES_P[:-1]

    for rank in ranks:
        sequences = read_sequences(gb, rank, args.source,
                                    args.species_from_file)
        for program in programs:
                alignment(sequences, program, rank, primers, args.num_threads)