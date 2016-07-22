library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

raxml_tree <- read.raxml("trees/RAxML_bipartitionsBranchLabels_esx4.newick")
ggtree(raxml_tree) + geom_label(aes(x=branch, label=bootstrap, fill=bootstrap)) + 
    geom_tiplab() +
    scale_fill_continuous(low='red', high='white') + 
    geom_treescale()