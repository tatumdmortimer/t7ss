library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

tree <- read.tree("trees/RAxML_bipartitionsBranchLabels.esx4_combine")
midpoint_rooted <- midpoint(tree)
reordered_tree <- reorder(midpoint_rooted)
p <- ggtree(midpoint_rooted) + geom_tiplab(size=3) + geom_treescale()
p

raxml_tree <- read.raxml("trees/RAxML_bipartitionsBranchLabels.esx4_midpointRoot.newick")
ggtree(raxml_tree) + geom_label(aes(x=branch, label=bootstrap, fill=bootstrap)) + 
    geom_tiplab() +
    scale_fill_continuous(low='red', high='white') + 
    geom_treescale()