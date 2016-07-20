# t7ss
Scripts and data associated with T7SS evolution paper

##aa_fastas
.faa files containing amino acids sequences for chromosomes and plasmids (output of prokka)

##alignments
Alignments used for phylogenetic analysis
####core_alignment.fasta
Concatenated alignment of amino acid sequences of all core genes of Actinobacteria used in this project
####eccA\_align\_trim.fasta
Alignment of eccA from plasmids
####eccB\_align\_trim.fasta
Alignment of eccB from plasmids
####espI\_align\_trim.fasta
Alignment of espI from plasmids
####esx4\_aa\_alignment.fasta
Concatenated alignment of amino acid sequences of ESX-4 genes from ESX-4 and ESX-4 bis loci in Actinobacteria
####esx\_aa\_alignmnet.fasta
Concatenated alignment of amino acid sequences of ESX genes for all ESX loci
####esx\_nucleotide\_alignment.fasta
Concatenated alignment of nucleotide sequences of ESX genes for all ESX loci
####plasmid\_ESX\_nucleotide\_align.fasta
Concatenated alignment of nucleotide sequences of ESX genes for plasmid-borne ESX loci
####tcpC\_align\_trim.fasta
Alignment of tcpC from plasmids
####ulcerans\_plasmid\_core\_alignment.fasta
Concatenated alignment of amino acid sequences of genes common to all _M. ulcerans_ plasmids.
####virB\_align\_trim.fasta
Alignment of virB from plasmids

##hyphy

####gblocks
Contains input and output of HyPhy analysis of alignment trimmed with GBlocks

####guidance
Contains input and output of HyPhy analysis of alignment trimmed with Guidance

####untrimmed
Contains input and output of HyPhy analysis of untrimmed alignment

##nucleotide_fastas
.ffn files containing nucleotide sequences for chromosomes and plasmids (output of prokka)

## orthomcl

####orthomcl_groups.txt
All orthologous groups output by OrthoMCL

####orthomcl\_core\_groups.txt
Orthologous groups corresponding to the core genome

####orthomcl\_esx\_groups.txt
Orthologous groups containing ESX genes

####orthomcl\_plasmids\_groups.txt
Orthologous groups in plasmids

## scripts

###checkPlasmidESX.py
Checks that all ESX genes identified on a plasmid are on the same component in the plasmidSPAdes assembly

Usage: checkPlasmidESX.py [-h] group plasmid

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
Downloads genomes from NCBI and names them genus\_species\_strain.fasta

Usage: getNCBIGenome.py [-h] accession email

####plotCoreGenomeTree.R
Creates Actinobacteria core phylogeny with ESX presence/absence matrix.

####plot_esx4_tree.R
Creates ESX-4 tree with bootstrap values colored according to support.

####run_prokka.sh
Runs prokka v 1.11 and compresses the output directory

Usage: run\_prokka.sh genus species strain genus\_species\_strain.fasta

####separateComponents.py
Uses orthomcl output to separate genes on the same component as ESX locus from other components and outputs new fasta files

Usage: separateComponents.py [-h] group plasmid

####sharedGeneContent.py
Calculates the number of shared genes between all pairs of plasmids in OrthoMCL output.

Usage: sharedGeneContent.py [-h] group categories

## submit_files
Submit files and DAGs to run scripts on HTCondor

#####get_genomes.dag

#####get_genomes.submit

#####prokka.dag

## tables

####esx\_gene\_names.txt
List of gene names of ESX genes in _M. tuberculosis_

####esx\_genes\_allCopies.txt
Table containing strain/locus names and locus tags for paralogs/orthologs of _eccA_, _eccB_, _eccC_, _eccD_, _eccE_, and _mycP_

####esx\_genes\_plasmids.txt
Same as above, but only plasmid-borne loci

####esxMatrix.txt
Presence/absence matrix of ESX loci in genomes in the core genome phylogeny

####finished\_plasmids\_accessions.txt
Accession numbers of finished plasmids

####genomeNames.txt
Table to match short names by OrthoMCL to full names

####group\_mtb\_genes.txt
OrthoMCL group names and the _M. tb_ ESX genes they contain

####mtb\_esx\_genes.txt
ESX genes in _M.tb_ with their Rv numbers and ESX locus

####plasmid_sra.txt
SRR/ERR accessions, genus, species, and strain information for assembled plasmids

####shared\_plasmid\_genes.txt
Comparison of number of shared genes between ESX containing plasmids and whetheror not they are in the same group (e.g. ESX-1P, ESX-2P, etc.)

####strains.txt
List of all genomes in analysis

####T7SS_accessions.tsv
Table with genomes and accessions used in analysis

##trees
####mrbayes
Contains trees from MrBayes phylogenetic analysis of eccA, eccB, virB, tcpC, and espI in plasmids
####RAxML_bipartitionsBranchLabels.core_longNames
Core genome phylogeny
####RAxML_bipartitionsBranchLabels_esx4.newick
ESX-4 phylogeny
####RAxML_bipartitionsBranchLabels.esx_combine_guidance
ESX phylogeny from finished genomes trimmed with Guidance
####RAxML_bipartitionsBranchLabels.esx_combine_trimmed
ESX phylogeney from finished genomes trimmed with Gblocks
####RAxML_bipartitionsBranchLabels.ulceransCore_combine
Phylogeny of _M. ulcerans_ plasmids



