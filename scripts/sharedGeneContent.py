#!/usr/bin/env python

import argparse
from collections import defaultdict

# Compares gene content between plasmids. Divides results by plasmids in the
# same group and plasmids from different groups

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Compare gene content')
    parser.add_argument("group", help="OrthoMCL groups file")
    parser.add_argument("categories", help="file describing plasmid groups")
    return parser.parse_args()


def get_orthomcl_genes(ortho_file):
    """Get gene content from OrthoMCL output"""
    geneContentDict = defaultdict(set)
    with open(ortho_file, 'r') as infile:
        for line in infile:
            line = line.strip().split(':')
            group = line[0]
            strains = [x.split('|')[0] for x in line[1].split()[1:]]
            for s in strains:
                geneContentDict[s].add(group)
    return geneContentDict

def get_categories(cat_file):
    """Get categories of plasmids from input file"""
    catDict = {}
    with open(cat_file, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            catDict[line[0]] = line[1]
    return catDict

def compare_gene_content(geneContentDict, catDict):
    strains = list(catDict.keys())
    with open("sharedGenes.txt", "w") as outfile:
        for i,strain1 in enumerate(strains):
            for strain2 in strains[i+1:]:
                sharedGenes = geneContentDict[strain1] & geneContentDict[strain2]
                sameGroup = "same"
                if catDict[strain1] != catDict[strain2]:
                    sameGroup = "different"
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    strain1, str(len(geneContentDict[strain1])), strain2,
                    str(len(geneContentDict[strain2])), str(len(sharedGenes)),
                    sameGroup))

args = get_args()
geneContent = get_orthomcl_genes(args.group)
catDict = get_categories(args.categories)
compare_gene_content(geneContent, catDict)
