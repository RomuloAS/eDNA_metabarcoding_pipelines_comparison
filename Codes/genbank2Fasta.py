#!/usr/bin/env python3
import os
import csv
import pathlib
import argparse
import datetime
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from Bio import Entrez
from ete3 import NCBITaxa

"""
Convert genbank format to fasta format to be used
in the pipelines execution.

Additionally to the conversion,
a taxonomy table with accession number + 
superkingdom,phylum,class,order,family,genus,species
is created when converting for Anacapa pipeline
and taxid table is created when converting
for SEQme pipeline.

"""

CHOICES_S = ["ncbi", "taxit", "taxdump"]

CHOICES_R = ["superkingdom", "phylum", "class",
            "order", "family", "genus",
            "species", "all"]

CHOICES_P = ["anacapa", "barque", "metabeat",
            "mifish", "seqme", "none", "all"]

PATH_MAIN = "{}_genbank2Fasta/".format(
                datetime.datetime.now().strftime("%d%m%Y_%H%M%S"))

LOG = "{}genbank2Fasta_error.log".format(PATH_MAIN)

PIPELINES_OUTPUT_FASTA = {"none": ">{species_}_{id}\n{seq}",
                        "anacapa": ">{id}\n{seq}",
                        "barque": ">{phylum}_{species_}\n{seq}",
                        "metabeat": ">{id}|{taxonid}|{species}\n{seq}",
                        "mifish": ">gb|{id}|{species_}\n{seq}",
                        "seqme": (">{id}\t{superkingdom};{kingdom};"
                                    "{phylum};{class};{order};"
                                    "{family};{genus};{species}\n{seq}")}

PIPELINES_OUTPUT_TAX = {"anacapa": ("{id}\t{superkingdom};{phylum};"
                                    "{class};{order};{family};"
                                    "{genus};{species}")}

LINEAGES = {}

TAXONOMIC_RANK = {"superkingdom": 0, "kingdom": 1, "phylum": 2,
                    "class": 3, "order": 4, "family": 5,
                    "genus": 6, "species": 7}
TAXONOMIC_HIERARCHY = {"kingdom": "superkingdom", "phylum": "kingdom",
                    "class": "phylum", "order": "class", "family": "order",
                    "genus": "family", "species": "genus"}

PATH_TO_TAXID = "{PATH_INFORMED}/TaxIDS.txt"
PATH_TO_TAXA = "{PATH_INFORMED}/Taxa.csv"
PATH_TO_DB = "{PATH_INFORMED}/ncbi_taxonomy.db"


def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    parser = argparse.ArgumentParser(
        description="Conversion from genbank to fasta format "\
                    "to be used in the execution of the pipeline(s).")
    parser.add_argument("gb", type=argparse.FileType("r"),
                   help="genbank file (.gb)")
    parser.add_argument('-sp', '--species_from_file', action='store_true',
                    help="Should it use species"\
                    " from file or download it from NCBI?")
    parser.add_argument("-p", "--pipeline", nargs="*", default=["none"],
                    choices=CHOICES_P, type = lambda s : s.lower(),
                    help="Pipeline convertion format, default: none")
    parser.add_argument("-r", "--rank", nargs="*", default=["superkingdom"],
                    choices=CHOICES_R, type = lambda s : s.lower(),
                    help="Taxonomic classification rank to be used"\
                    " to separate the groups, default: superkingdom")
    parser.add_argument("-s", "--source", nargs="?", const="taxdump",
                    default="taxdump", choices=CHOICES_S,
                    type = lambda s : s.lower(),
                    help="Source to be used to collect"\
                    " the info about the taxonomic rank, default: taxdump")
    parser.add_argument('-t', '--path_to_taxid_files', nargs='?', type = str,
                    const=PATH_MAIN, default=PATH_MAIN,
                    help="Path to taxit files.")
    
    return parser.parse_args()

def parse_taxID(gb):
    """Parse taxon ids from genbank file

    This function uses Biopython library to parse
    taxon ids and create a file with them.

    Parameters
    -------------
    gb: io.TextIOWrapper
        Genbank file

    """

    tax_ids = set()

    for r in tqdm(gb, desc="Reading sequences"):
        record = gb.get(r)
        taxonid = record.features[0].qualifiers["db_xref"][0].split(":")[1]
        tax_ids.add(taxonid)

    with open(PATH_TO_TAXID, "w") as out_taxids:
        out_taxids.write("\n".join(tax_ids))

def taxit():
    """Download a create taxonomic database using taxit

    This function executes taxit to download taxonomy database
    and creates a table with the taxonomic lineages.
    """

    print("Downloading taxit database ...")
    subprocess.call("taxit new_database {PATH_TO_DB} -p {PATH_DOWNLOAD}".format(
                    PATH_TO_DB=PATH_TO_DB, 
                    PATH_DOWNLOAD=PATH_MAIN) , shell=True)

    print("Creating tax table ...")
    subprocess.call("taxit taxtable  {PATH_TO_DB}".format(PATH_TO_DB=PATH_TO_DB) + 
                    " -f {PATH_TO_TAXID}".format(PATH_TO_TAXID=PATH_TO_TAXID) + 
                    " -o {PATH_TO_TAXA}".format(PATH_TO_TAXA=PATH_TO_TAXA),
                    shell=True)

def parse_taxa():
    """Parse taxonomic information from Taxa.csv

    This function opens Taxa.csv file to parse
    tax id and lineage rank.
    
    Returns
    -------------
    tax_rank_id: dict
        Taxonomic rank and id
    """

    with open(PATH_TO_TAXA) as taxa_input:
        taxa = csv.DictReader(taxa_input)
        tax_rank_id = {row["tax_id"]:row for row in taxa}

    return tax_rank_id

def get_tax_lineage(taxonid, source, tax_rank_id={}):
    """Return taxonomy lineage information

    This function uses either Biopython library to connect
    NCBI database and search for taxonomy information
    or searches the information locally by using ete3 taxdump
    file or taxit program to create sql version of it.

    Parameters
    -------------
    taxonid : string
        Taxonomic id of the species
    source : string
        Source to be used to collect the info about the taxonid
    tax_rank_id: dict
        Taxonomic rank and id

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

        if source == "taxit":
            lineage = {level:tax_rank_id[tax_rank_id[
                        taxonid][level]]["tax_name"]
                        for level in CHOICES_R[:-1]}

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

def read_sequences(gb, pipeline, rank, source, species_from_file, tax_rank_id={}):
    """Read the genbank file and parse the sequences
    based on the taxonomic rank

    This function uses Biopython library to scan the genbank file
    and parse the sequences based on the taxonomic rank.

    Parameters
    -------------
    gb : io.TextIOWrapper
        A genbank file
    pipeline: string
        The pipeline format of the fasta format
    rank: string
        Taxonomic rank
    source : string
        Source to be used to collect the info about the taxonid
    species_from_file: bool
        Indicate if species name from file should be used
    tax_rank_id: dict
        Taxonomic rank and id

    Returns
    -------------
    sequences: dictionary
        A dictionary with key representing ranks and values
        representing sequences in fasta format with species name
        as the header of the sequence
    tax_tables: dictionary
        A dictionary with key representing ranks and values
        representing taxonomic tables if the pipeline
        is either anacapa or seqme 

    """

    sequences = {}
    tax_tables = {}
    sequence_info = {}
    taxid_table = {}

    for r in tqdm(gb, desc="Reading sequences"):
        record = gb.get(r)
        tax = []

        sequence_info["id"] = record.id
        sequence_info["seq"] = record.seq

        for feature in record.features:
            if feature.type == "source" and \
                "taxon" in feature.qualifiers["db_xref"][0]:
                taxonid = feature.qualifiers["db_xref"][0].split(":")[1]
                sequence_info["taxonid"] = taxonid

        lineage = get_tax_lineage(taxonid, source, tax_rank_id)

        if species_from_file:
            lineage["species"] = record.features[0].qualifiers[
                                        "organism"][0].lstrip().rstrip()
        else:
            if "species" not in lineage:
                lineage["species"] = record.features[0].qualifiers[
                                        "organism"][0].lstrip().rstrip()

        lineage["species_"] = lineage["species"].replace(" ", "_")
        sequence_info.update(lineage)

        if pipeline in "anacapa":
                tax = PIPELINES_OUTPUT_TAX[pipeline].format(**sequence_info)
        elif pipeline == "seqme":

            for lin in lineage:
                if lin not in TAXONOMIC_RANK:
                    continue

                if lineage[rank] not in taxid_table:
                    taxid_table[lineage[rank]] = [
                            {"Eukaryota":"0*Eukaryota*-1*0*superkingdom"},
                            {"Eukaryota": 0}, 1
                        ]
                    tax.append(taxid_table[lineage[rank]][0]["Eukaryota"])

                if lineage[lin] not in taxid_table[lineage[rank]][0]:
                    taxid_table[lineage[rank]][1][
                            lineage[lin]] = taxid_table[lineage[rank]][2]
                    taxid_table[lineage[rank]][0][
                            lineage[lin]] = "{number}*{tax}*{b_tax}"\
                                    "*{taxid}*{lineage}".format(
                                    number=taxid_table[lineage[rank]][2], 
                                    tax=lineage[lin],
                                    b_tax=taxid_table[lineage[rank]][1][lineage[
                                            TAXONOMIC_HIERARCHY[lin]]],
                                    taxid=TAXONOMIC_RANK[lin],
                                    lineage=lin)

                    taxid_table[lineage[rank]][2] += 1
                    tax.append(taxid_table[lineage[rank]][0][lineage[lin]])

            tax = "\n".join(tax)

        try:
            sequence = PIPELINES_OUTPUT_FASTA[pipeline].format(**sequence_info)

            if lineage[rank] not in sequences:
                sequences[lineage[rank]] = [sequence]
                if tax:
                    tax_tables[lineage[rank]] = [tax]
            else:
                sequences[lineage[rank]].append(sequence)
                if tax:
                    tax_tables[lineage[rank]].append(tax)

        except:
            with open(LOG, "a") as log:
                print("\nRank '{}' not found for organism '{}', taxonid '{}'".format(
                    rank, lineage["species"], taxonid), file=log)

    return sequences, tax_tables

def save_fasta(sequences, pipeline, rank):
    """Save sequences to file

    This function saves each group of sequence
    in the dictionary to fasta file format
    based on rank grouping.

    sequences : dictionary
        A dictionary with sequences
    pipeline: string
        The pipeline format of the fasta format
    rank: string
        Taxonomic rank

    """

    path_fasta = "{}{}/{}".format(PATH_MAIN, pipeline, rank)
    pathlib.Path(path_fasta).mkdir(parents=True, exist_ok=True)

    for seq in tqdm(sequences, desc="Saving FASTA"):
        with open("{PATH_FASTA}/{FASTA_NAME}.fasta".format(
                                    PATH_FASTA=path_fasta,
                                    FASTA_NAME=seq.replace(" ", "_")
                                    ), "w") as fasta_file:
            fasta_file.write("\n".join(sequences[seq]))

def save_tax_tables(tax_tables, pipeline, rank):
    """Save tax tables to file
    
    This function saves each group of tax table
    in the dictionary to a text file
    based on rank grouping. Tax tables are
    only created if the pipeline variable
    is equal to anacapa or seqme.

    tax_tables : dictionary
        A dictionary with tax tables
    pipeline: string
        The pipeline format of the fasta format
    rank: string
        Taxonomic rank

    """

    path_tax_tables = "{}{}/{}".format(PATH_MAIN, pipeline, rank)
    pathlib.Path(path_tax_tables).mkdir(parents=True, exist_ok=True)

    for tax_table in tqdm(tax_tables, desc="Saving tax table"):
        with open("{PATH_TAX_TABLE}/{TAX_TABLE_NAME}.txt".format(
                                    PATH_TAX_TABLE=path_tax_tables,
                                    TAX_TABLE_NAME=tax_table.replace(
                                                    " ", "_")
                                    ), "w") as tax_table_file:
            tax_table_file.write("\n".join(tax_tables[tax_table]))


if __name__ == "__main__":
    args = getArguments()
    
    pathlib.Path(PATH_MAIN).mkdir(parents=True, exist_ok=True)
    gb = SeqIO.index(args.gb.name, "genbank")
    with open(LOG, "w"): pass

    tax_rank_id = {}
    if args.source == "taxit":

        PATH_TO_TAXID = "{PATH_INFORMED}/TaxIDS.txt".format(
                        PATH_INFORMED=args.path_to_taxid_files)
        PATH_TO_TAXA = "{PATH_INFORMED}/Taxa.csv".format(
                        PATH_INFORMED=args.path_to_taxid_files)
        PATH_TO_DB = "{PATH_INFORMED}/ncbi_taxonomy.db".format(
                        PATH_INFORMED=args.path_to_taxid_files)

        if not os.path.isfile(PATH_TO_TAXA):
            parse_taxID(gb)
            taxit()

        tax_rank_id = parse_taxa()

    ranks = args.rank
    if "all" in args.rank:
        ranks = CHOICES_R[:-1]

    pipelines = args.pipeline
    if "all" in args.pipeline:
        pipelines = CHOICES_P[:-1]

    for pipeline in pipelines:
        print("Pipeline: {}".format(pipeline))
        for rank in ranks:
            print("Rank: {}".format(rank))

            sequences, tax_tables = read_sequences(gb, pipeline, 
                                                rank, args.source,
                                                args.species_from_file,
                                                tax_rank_id)
            
            save_fasta(sequences, pipeline, rank)
            if tax_tables:
                save_tax_tables(tax_tables, pipeline, rank)