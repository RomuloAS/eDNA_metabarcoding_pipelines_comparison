#!/usr/bin/env python3
import os
import csv
import pathlib
import argparse
import datetime
import subprocess
from Bio import SeqIO
from Bio import Entrez
from tqdm import tqdm
from ete3 import NCBITaxa

"""
The SATIVA Algorithm is used for taxonomically
mislabelled sequences identification and to
suggest corrections.
"""

CHOICES_S = ["ncbi", "taxit", "taxdump"]

TAX_LEVELS = ["superkingdom","phylum","class","order","family","genus","species"]
LINEAGES = {}

PATH_MAIN = "{}_sativa/".format(datetime.datetime.now().strftime("%d%m%Y_%H%M%S"))
PATH_TO_SATIVA_TAX = "{PATH_MAIN}Sativa.tax".format(PATH_MAIN=PATH_MAIN) 
PATH_TO_SATIVA = "sativa/"

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
        description="Identification of taxonomically"\
                    " mislabelled sequences")
    parser.add_argument("gb", type=argparse.FileType("r"),
                   help="genbank format file (.gb)")
    parser.add_argument("phy", type=argparse.FileType("r"),
                   help="PHYLIP multiple sequence alignment format file (.phy)")
    parser.add_argument("-s", "--source", nargs="?", const="taxdump", default="taxdump",
                    choices=CHOICES_S, type = lambda s : s.lower(),
                    help="Source to be used to collect"\
                    " the info about the taxonomic rank, default: taxdump")
    parser.add_argument('-p', '--path_to_sativa', nargs='?', type = str,
                    const=PATH_TO_SATIVA, default=PATH_TO_SATIVA,
                    help="Path to sativa code.")
    parser.add_argument('-t', '--path_to_taxid_files', nargs='?', type = str,
                    const=PATH_MAIN, default=PATH_MAIN,
                    help="Path to taxit files.")
    
    return parser.parse_args()

def download_and_install_sativa():
    """Check if sativa is installed

    This function checks if the path to sativa exists
    and if sativa was installed.

    """
        
    print("Downloading sativa ...")
    subprocess.call("git clone --recursive" \
                    " https://github.com/"\
                    "amkozlov/sativa.git {path_to_sativa}".format(
                    path_to_sativa=PATH_TO_SATIVA),
                    shell=True)

    print("Installing sativa ...")
    subprocess.call("bash {path_to_sativa}/install.sh".format(
                    path_to_sativa=PATH_TO_SATIVA),
                    shell=True)

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
                        for level in TAX_LEVELS}

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
    gb = SeqIO.index(gb.name, "genbank")

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

    print("Downloading database ...")
    subprocess.call("taxit new_database {PATH_TO_DB} -p {PATH_DOWNLOAD}".format(
                    PATH_TO_DB=PATH_TO_DB, 
                    PATH_DOWNLOAD=PATH_MAIN) , shell=True)

    print("Creating tax table ...")
    subprocess.call("taxit taxtable  {PATH_TO_DB}".format(PATH_TO_DB=PATH_TO_DB) + 
                    " -f {PATH_TO_TAXID}".format(PATH_TO_TAXID=PATH_TO_TAXID) + 
                    " -o {PATH_TO_TAXA}".format(PATH_TO_TAXA=PATH_TO_TAXA),
                    shell=True)

    print("DONE!")

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

def tax4Sativa(gb, source, tax_rank_id={}):
    """Generate a taxonomic file for sativa
    using NCBI, taxa.csv from taxit, or taxdump with ete3

    This function creates a file with the taxon id,
    and the tax levels for ech taxon id.
    
    Parameters
    -------------
    gb: io.TextIOWrapper
        Genbank file
    source : string
        Source to be used to collect the info about the taxonid
    tax_rank_id: dict
        Taxonomic rank and id
    """

    sativa_taxes = []

    gb = SeqIO.index(gb.name, "genbank")

    for r in tqdm(gb, desc="Creating Sativa.tax"):
        record = gb.get(r)
        taxonid = record.features[0].qualifiers["db_xref"][0].split(":")[1]
        lineage = get_tax_lineage(taxonid, source, tax_rank_id)

        sativa_tax = "{}\t".format(record.id)
        for level in TAX_LEVELS:
            if level not in lineage:
                sativa_tax += "unknown;"
                continue

            sativa_tax += "{};".format(lineage[level])

        sativa_taxes.append(sativa_tax[:-1])

    with open("{PATH_TO_SATIVA_TAX}".format(
                PATH_TO_SATIVA_TAX=PATH_TO_SATIVA_TAX
                ), "w") as sativa_tax_output:
        sativa_tax_output.write("\n".join(sativa_taxes))

def sativa(phy):
    """Run sativa for identification of taxonomically
    mislabelled sequences

    This function executes SATIVA algorithm to identify
    taxonomically mislabelled sequences.

    phy: io.TextIOWrapper
        PHYLIP multiple sequence alignment format
    """
    
    print("SATIVA ...")
    subprocess.call("python {PATH_TO_SATIVA}/sativa.py".format(
                    PATH_TO_SATIVA=PATH_TO_SATIVA) +
                    " -s {PATH_TO_PHY} -t {PATH_TO_SATIVA_TAX}".format(
                    PATH_TO_PHY=phy.name, PATH_TO_SATIVA_TAX=PATH_TO_SATIVA_TAX) +
                    " -x zoo -n 12S -o {PATH_MAIN}sativa_result/".format(
                    PATH_MAIN=PATH_MAIN) + 
                    " -T 10 -v" , shell=True)
    print("DONE!")
            

if __name__ == "__main__":
    args = getArguments()

    pathlib.Path("{PATH_MAIN}sativa_result/".format(
                PATH_MAIN=PATH_MAIN)).mkdir(
                parents=True, exist_ok=True)
    PATH_TO_SATIVA = args.path_to_sativa

    if not os.path.isfile("{path_to_sativa}/sativa.py".format(
                            path_to_sativa=PATH_TO_SATIVA)):
        download_and_install_sativa()

    if args.source == "ncbi":

        tax4Sativa(args.gb, args.source)
        sativa(args.phy)

    elif args.source == "taxit":

        PATH_TO_TAXID = "{PATH_INFORMED}/TaxIDS.txt".format(
                        PATH_INFORMED=args.path_to_taxid_files)
        PATH_TO_TAXA = "{PATH_INFORMED}/Taxa.csv".format(
                        PATH_INFORMED=args.path_to_taxid_files)
        PATH_TO_DB = "{PATH_INFORMED}/ncbi_taxonomy.db".format(
                        PATH_INFORMED=args.path_to_taxid_files)

        if not os.path.isfile(PATH_TO_TAXA):
            parse_taxID(args.gb)
            taxit()

        tax_rank_id = parse_taxa()
        tax4Sativa(args.gb, args.source, tax_rank_id)
        sativa(args.phy)

    elif args.source == "taxdump":
        
        tax4Sativa(args.gb, args.source)
        sativa(args.phy)