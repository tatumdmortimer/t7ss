library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

tree <- read.tree("RAxML_bestTree.core_bestTree")
midpoint_rooted <- midpoint(tree)
reordered_tree <- reorder(midpoint_rooted)
p <- ggtree(reordered_tree) + geom_tiplab(size=2)
dd <- read.table("esxMatrix_shortNames.txt", sep="\t", stringsAsFactor=F)
dd[dd=="no"] <- ""
p <- gheatmap(p, dd, offset = 0.5, width = 0.5, colnames_position="top", font.size=1) + 
    theme(legend.position="none") +
    scale_fill_manual(values=c("black"))
p
## can't find a way to rotate the labels, so may have to be done in illustrator
