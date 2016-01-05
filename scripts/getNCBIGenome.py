#!/usr/bin/env python

import argparse
from Bio import Entrez
from Bio import SeqIO


# This script downloads genomes from NCBI using accession number

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Download genome')
    parser.add_argument("accession", help="NCBI Accession Number")
    parser.add_argument("email", help="E-mail address")
    return parser.parse_args()

args = get_args()

Entrez.email = args.email
handle = Entrez.efetch(db="nucleotide", id=args.accession, 
    rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

name = record.description.split(",")[0].split(" ")[1:]
genus = name[0]
species = name[1]
if name[2] == "strain" or name[2] == "str." or name[2] == "subsp.":
    strain = name[3].replace('/', '')
elif name[2] == "ATCC":
    strain = name[2] + name[3].replace('/', '')
else:
    strain = name[2].replace('/', '')
print genus, species, strain

SeqIO.write(record, "_".join([genus,species,strain]) + ".fasta", "fasta")
