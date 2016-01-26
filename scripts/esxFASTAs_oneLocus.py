#!/usr/bin/python

import argparse
import os
from Bio import SeqIO

################################################################################
# This script creates fasta files of unaligned orthologs and paralogs from   
# input files describing gene names and .faa and .fna files of annotated genomes
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
    parser = argparse.ArgumentParser(description='Create fasta of esx genes')
    parser.add_argument("esx", help="ESX genes file", action=FullPaths,
        type=is_file)
    parser.add_argument("nucleotide", 
        help="Directory with .fasta files of gene sequences in nucleotides",
        action=FullPaths, type=is_dir)
    parser.add_argument("protein", 
        help="Directory with .fasta files of protein sequences in amino acids",
        action=FullPaths, type=is_dir)
    return parser.parse_args()

def concatenate_sequences(seq1, seq2, newName):
    seq1.seq = seq1.seq + seq2.seq
    seq1.id = newName
    return seq1

def make_unaligned_fasta(dnaDirectory, protDirectory, geneList, geneName):
    """ Reads through files in provided directory to find gene sequences that
    match the proteins in input files for each esx locus"""
    species = ["_".join(x.split("_")[:-1]) for x in geneList]
    dnaSeqDict = {}
    protSeqDict = {}
    for s in species:
        with open(dnaDirectory + '/' + s + ".ffn", "r") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                dnaSeqDict[record.id] = record
        with open(protDirectory + '/' + s + ".faa", "r") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                protSeqDict[record.id] = record
    dna_out = open(geneName + '.ffn', 'w')
    prot_out = open(geneName + '.faa', 'w')
    dna_records = []
    prot_records = []
    for gene in geneList:
        if '/' in gene:
            gene1 = gene.split('/')[0]
            gene2 = gene.split('/')[0][:-2] + gene.split('/')[1]
            dna_records.append(concatenate_sequences(dnaSeqDict[gene1],
                dnaSeqDict[gene2], gene))
            prot_records.append(concatenate_sequences(protSeqDict[gene1],
                protSeqDict[gene2], gene))
        else:
            dna_records.append(dnaSeqDict[gene])
            prot_records.append(protSeqDict[gene])
    SeqIO.write(dna_records, dna_out, 'fasta')
    SeqIO.write(prot_records, prot_out, 'fasta')

args = get_args()

eccA = []
eccB = []
eccC = []
eccD = []
eccE = []
mycP = []

with open(args.esx, 'r') as infile:
    for i, line in enumerate(infile):
        line = line.strip().split()
        eccA.append(line[1])
        eccB.append(line[2])
        eccC.append(line[3])
        eccD.append(line[4])
        eccE.append(line[5])
        mycP.append(line[6])

eccA = [x for x in eccA if x != '-']
eccB = [x for x in eccB if x != '-']
eccC = [x for x in eccC if x != '-']
eccD = [x for x in eccD if x != '-']
eccE = [x for x in eccE if x != '-']
mycP = [x for x in mycP if x != '-']

for j,k in [(eccA, 'eccA'), (eccB, 'eccB'), (eccC, 'eccC'), (eccD, 'eccD'),
    (eccE, 'eccE'), (mycP, 'mycP')]:
    make_unaligned_fasta(args.nucleotide, args.protein, j, k)
