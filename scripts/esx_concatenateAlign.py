#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

################################################################################
# This script creates a concatenated alignment from aligned esx genes
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Create concatenated alignment\
 of ESX genes')
    parser.add_argument("eccA", help="eccA alignment", action=FullPaths,
        type=is_file)
    parser.add_argument("eccB", help="eccB alignment", action=FullPaths,
        type=is_file)
    parser.add_argument("eccC", help="eccC alignment", action=FullPaths,
        type=is_file)
    parser.add_argument("eccD", help="eccD alignment", action=FullPaths,
        type=is_file)
    parser.add_argument("eccE", help="eccE alignment", action=FullPaths,
        type=is_file)
    parser.add_argument("mycP", help="mycP alignmnet", action=FullPaths,
        type=is_file)
    parser.add_argument("description", help="file with loci and genes names",
        action=FullPaths, type=is_file)
    return parser.parse_args()

def concatenate_sequences(seqs, newName):
    seqs[0].seq = seqs[0].seq + seqs[1].seq + seqs[2].seq + seqs[3].seq + seqs[4].seq + seqs[5].seq
    seqs[0].id = newName
    return seqs[0]

def check_align_get_length(seq_dict):
    length = 0
    for rec in seq_dict:
        if length == 0:
            length = len(seq_dict[rec].seq)
        elif len(seq_dict[rec].seq) != length:
            print("At least one alignment has records of different lengths")
            sys.exit()
    return length


args = get_args()

eccA_dict = SeqIO.to_dict(SeqIO.parse(open(args.eccA, "r"), "fasta"))
eccB_dict = SeqIO.to_dict(SeqIO.parse(open(args.eccB, "r"), "fasta"))
eccC_dict = SeqIO.to_dict(SeqIO.parse(open(args.eccC, "r"), "fasta"))
eccD_dict = SeqIO.to_dict(SeqIO.parse(open(args.eccD, "r"), "fasta"))
eccE_dict = SeqIO.to_dict(SeqIO.parse(open(args.eccE, "r"), "fasta"))
mycP_dict = SeqIO.to_dict(SeqIO.parse(open(args.mycP, "r"), "fasta"))

lengths = [check_align_get_length(x) for x in [eccA_dict, eccB_dict, eccC_dict,
    eccD_dict, eccE_dict, mycP_dict]]

concatenated_seqs = []

with open(args.description, "r") as infile:
    for line in infile:
        line = line.strip().split()
        locus_name = line[0]
        gene_names = line[1:]
        genes = []
        for i,g in enumerate(gene_names):
            if g == "-":
                genes.append("-")
            elif i == 0:
                genes.append(eccA_dict[g])
            elif i == 1:
                genes.append(eccB_dict[g])
            elif i == 2:
                genes.append(eccC_dict[g])
            elif i == 3:
                genes.append(eccD_dict[g])
            elif i == 4:
                genes.append(eccE_dict[g])
            elif i == 5:
                genes.append(mycP_dict[g])
        gaps_index = [i for i,j in enumerate(genes) if j == "-"]
        for gap in gaps_index:
            new_record = SeqRecord(Seq("".join(["-"]*lengths[gap])),
            id="placeholder")
            genes[gap] = new_record
        concatenated_seqs.append(concatenate_sequences(genes, locus_name))

SeqIO.write(concatenated_seqs, "ESX_concatenated_alignment.fasta", "fasta")




