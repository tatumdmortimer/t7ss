#!/usr/bin/env python

import argparse
from collections import defaultdict
from Bio import SeqIO


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


def get_component(geneDict, plasmidFile):
    """Check that all genes are on the same component and output new annotation
    files for that component only"""
    componentDict = {}
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
            componentDict[strain] = components
    return componentDict

def create_fastas(componentDict, plasmidFile):
    """Create new nucleotide and amino acid fasta files that are separate
    for each component from the assembly with ESX genes"""
    with open(plasmidFile) as infile:
        for line in infile:
            line = line.strip().split()
            strain = line[0]
            aa = 'newAnnotations/{0}_plasmid/{0}_plasmid.faa'.format(strain)
            nuc = 'newAnnotations/{0}_plasmid/{0}_plasmid.ffn'.format(strain)
            gff = 'newAnnotations/{0}_plasmid/{0}_plasmid.gff'.format(strain)
            aa_dict = SeqIO.to_dict(SeqIO.parse(aa, "fasta"))
            nuc_dict = SeqIO.to_dict(SeqIO.parse(nuc, "fasta"))
            for c in componentDict[strain]:
                new_aa = []
                new_nuc = []
                with open(gff, 'r') as gff_file:
                    for line in gff_file:
                        if line[0] != '#':
                            if line[0] == '>':
                                break
                            line = line.strip().split()
                            component = line[0][-2]
                            gene = line[8].split(';')[0][3:]
                            if component == c:
                                try:
                                    new_aa.append(aa_dict[gene])
                                    new_nuc.append(nuc_dict[gene])
                                except KeyError:
                                    print("{0} does not exist".format(gene))
                SeqIO.write(new_aa, "{0}_{1}.faa".format(strain, c), "fasta")
                SeqIO.write(new_nuc, "{0}_{1}.ffn".format(strain, c), "fasta")

args = get_args()
geneDict = get_orthomcl_genes(args.group)
componentDict = get_component(geneDict, args.plasmid)
create_fastas(componentDict, args.plasmid)
