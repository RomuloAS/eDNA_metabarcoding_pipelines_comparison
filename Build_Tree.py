#!/usr/bin/env python3
import os
import argparse
import datetime
import subprocess
from ete3 import PhyloTree, TreeStyle

"""
Build Phylogenetic Tree using raxmlHPC-PTHREADS-SSE3

From the alignment in fasta format and tree in newick format,
a pdf with the tree and alignment (side by side) will be generated.
"""

RAXML = ("raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA" 
        " -n {output_file} -p 765 -s {input_file}"
        " -T 10 -x 498 -N 100")

CHOICES_F = ["pdf", "svg", "png", "jpg"]

def getArguments():
    """Get arguments from terminal

    This function gets arguments from terminal via argparse

    Returns
    -------------
    arguments: Namespace
        Namespace object with all arguments
    """

    parser = argparse.ArgumentParser(
        description="Build the Phylogenetic Tree and save it in a file.")
    parser.add_argument("aln", type=argparse.FileType("r"),
                   help="alignment file")
    parser.add_argument("-f", "--format", nargs="?", const="jpg", default="jpg",
                    choices=CHOICES_F, type = lambda s : s.lower(),
                    help="Format to save the phylogenetic tree, default: jpg")
    parser.add_argument('-s', '--show', action='store_true',
                    help='Show ETE tree Browser')
    
    return parser.parse_args()

def raxml(aln, basename):
    """Build a phylogenetic tree based on alignment provided
    using raxmlHPC-PTHREADS-SSE3

    This function executes RAxML program the create a
    phylogenetic tree using GTRGAMMA substitution model 

    Parameters
    -------------
    aln: io.TextIOWrapper
        Alignment file
    basename: string
        Basename of the original alignment file
    """

    print("Building tree ...")
    output_file = "{BASENAME}.raxml".format(BASENAME=basename)

    subprocess.call(RAXML.format(output_file=output_file,
                                input_file=aln.name),
                                shell=True)
    print("RAxML DONE!")

def build_tree(aln, tree, basename, show, output_format):
    """Build phylogenetic tree from files

    This function creates a file with the phylogenetic tree and alignment
    from the fasta multiple alignment file and the tree in newick format.

    Parameters
    -------------
    aln: string
        Alignment string in fasta format
    tree: string
        Tree string in newick format
    basename: string
        Basename of the original alignment file
    show: boolean
        Show ETE tree browser (yes/no)
    output_format: string
        Format of the output
    """

    if tree[-1] != ";":
        genetree = PhyloTree("{};".format(tree))
    else:
        genetree = PhyloTree(tree)

    ts = TreeStyle()
    ts.show_leaf_name = False

    new_tree = "{BASENAME}_Tree.{FORMAT}".format(
                                BASENAME=basename,
                                FORMAT=output_format)
    new_tree_aln = "{BASENAME}_Tree_aln.{FORMAT}".format(
                                BASENAME=basename,
                                FORMAT=output_format)
    
    if show:        
        genetree.render(new_tree, tree_style=ts)
        genetree.link_to_alignment(aln)
        genetree.render(new_tree_aln, tree_style=ts)
        genetree.show(tree_style=ts)
    else:
        genetree.render(new_tree, tree_style=ts)
        genetree.link_to_alignment(aln)
        genetree.render(new_tree_aln, tree_style=ts)

if __name__ == "__main__":
    args = getArguments()

    basename = os.path.splitext(os.path.basename(
                                args.aln.name))[0]

    raxml(args.aln, basename)

    tree_file = "RAxML_bestTree.{BASENAME}.raxml".format(
                                    BASENAME=basename)

    aln = args.aln.read()
    with open(tree_file) as tree_file:
        tree = tree_file.read()

    build_tree(aln, tree, basename, args.show, args.format)