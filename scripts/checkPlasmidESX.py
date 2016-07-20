#!/usr/bin/env python

import argparse
from collections import defaultdict


# This gets conserved gene names from orthomcl output, checks that those genes
# are on the same component from plasmidSPAdes assembly, and outputs a file
# listing gene names.

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Check plasmid ESX')
    parser.add_argument("group", help="OrthoMCL groups file")
    parser.add_argument("plasmid", help="File containing plasmid info")
    return parser.parse_args()


def get_orthomcl_genes(ortho_file):
    """Get gene names for conserved ESX genes from OrthoMCL output"""
    eccA = ['core1161']
    eccB = ['core1077']
    eccC = ['core1132']
    eccD = ['core1266', 'core2180']
    eccE = ['core1198']
    mycP = ['core1159']
    groupDict = {}
    with open(ortho_file, 'r') as infile:
        for line in infile:
            line = line.strip().split(':')
            groupDict[line[0]] = line[1].split()[1:]
    geneDict = defaultdict(lambda: ['-', '-', '-', '-', '-', '-'])
    for g in groupDict:
        genes = groupDict[g]
        for gene in genes:
            strain = gene.split('|')[1][:-6]
            if g in eccA:
                geneDict[strain][0] = gene.split('|')[1]
            elif g in eccB:
                geneDict[strain][1] = gene.split('|')[1]
            elif g in eccC:
                geneDict[strain][2] = gene.split('|')[1]
            elif g in eccD:
                geneDict[strain][3] = gene.split('|')[1]
            elif g in eccE:
                geneDict[strain][4] = gene.split('|')[1]
            elif g in mycP:
                geneDict[strain][5] = gene.split('|')[1]
    return geneDict


def check_component(geneDict, plasmidFile):
    """Check that all genes are on the same component and output new annotation
    files for that component only"""
    with open(plasmidFile, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            strain = line[0]
            esx_genes = geneDict[strain]
            components = set()
            gff = 'newAnnotations/{0}_plasmid/{0}_plasmid.gff'.format(strain)
            with open(gff, 'r') as gff_file:
                for line in gff_file:
                    if line[0] != '#':
                        if line[0] == '>':
                            break
                        line = line.strip().split()
                        component = line[0][-2]
                        gene = line[8].split(';')[0][3:]
                        if gene in esx_genes:
                            components.add(component)
            if len(components) > 1:
                print('WARNING: ESX genes for strain {0} on components {1}!'.format(strain, ', '.join(components)))


args = get_args()
geneDict = get_orthomcl_genes(args.group)
check_component(geneDict, args.plasmid)
