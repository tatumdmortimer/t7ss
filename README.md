# t7ss
Scripts and data associated with T7SS evolution paper

## scripts

####concatenateEccC1.py
Concatenates fasta format sequences of eccCa1 and eccCb1
Usage: concatenateEccC1.py [-h] eccCa1 eccCb1

####esxFASTAs_oneLocus.py
Creates fasta files of nucleotide and amino acid sequences for alignment of ESX genes. Requires input file formatted as follows:
```
species1_locus1 eccA_gene_name  eccB_gene_name  eccC_gene_name  eccD_gene_name  eccE_gene_name  mycP_gene_name
species2_locus1 eccA_gene_name  eccB_gene_name  eccC_gene_name  eccD_gene_name  eccE_gene_name  mycP_gene_name
species2_locus2 eccA_gene_name  eccB_gene_name  eccC_gene_name  eccD_gene_name  eccE_gene_name  mycP_gene_name
```
Usage: esxFASTAs_oneLocus.py [-h] esx nucleotide protein

####esxFASTAs.py
Creates fasta files of nucleotide and amino acid sequences for alignmnent. Creates separate files for each ESX locus. Each ESX locus requires an input file as described for esxFASTAs_oneLocus.py.

Usage: esxFASTAs.py [-h] esx1 esx2 esx3 esx4 esx5 nucleotide protein

####esx_treescape.R
Performs treescape analysis on MrBayes output

####getNCBIGenome.py
Downloads genomes from NCBI and names them {genus}_{species}_{strain}.fasta
Usage: getNCBIGenome.py [-h] accession email

####plotCoreGenomeTree.R
Creates Actinobacteria core phylogeny with ESX presence/absence matrix.

####plot_esx4_tree.R
Creates ESX-4 tree with bootstrap values colored according to support.

####run_prokka.sh
Runs prokka v 1.11 and compresses the output directory

Usage: run_prokka.sh genus species strain genus_species_strain.fasta

## submit_files

#####get_genomes.dag

#####get_genomes.submit

#####prokka.dag

## tables

#####esx_gene_names.txt

#####group_mtb_genes.txt

#####mtb_esx_genes.txt

#####strains.txt

#####T7SS_accessions.tsv

#####shared_plasmid_genes.txt
Comparison of number of shared genes between ESX containing plasmids and whetheror not they are in the same group (e.g. ESX-1P, ESX-2P, etc.)

## orthomcl

#####orthomcl_groups.txt

#####orthomcl_core_groups.txt

#####orthomcl_esx_groups.txt
