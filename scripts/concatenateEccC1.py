#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# This script concatenates the eccC1a and eccC1b sequences

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Concatenate eccC1')
    parser.add_argument("eccCa1", help="fasta file of eccCa1 sequences")
    parser.add_argument("eccCb1", help="fasta file of eccCb1 sequences")
    return parser.parse_args()

args = get_args()

a_seqs = {}
b_seqs = {}

for a_rec in SeqIO.parse(open(args.eccCa1), "fasta"):
    strain = "_".join(a_rec.id.split("_")[:-1])
    a_seqs[strain] = str(a_rec.seq)

for b_rec in SeqIO.parse(open(args.eccCb1), "fasta"):
    strain = "_".join(b_rec.id.split("_")[:-1])
    b_seqs[strain] = str(b_rec.seq)

a_set = set(a_seqs.keys())
b_set = set(b_seqs.keys())

print len(a_set - b_set)
print len(b_set - a_set)
strains = list(a_set & b_set)

cat_seqs = []

for s in strains:
    record = SeqRecord(Seq(a_seqs[s] + b_seqs[s]), id=s,
                       description="concatenated eccC1")
    cat_seqs.append(record)

SeqIO.write(cat_seqs, "concatenatedEccC1.fasta", "fasta")
