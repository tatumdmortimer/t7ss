library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

tree <- read.tree("trees/RAxML_bipartitionsBranchLabels.core_longNames")
midpoint_rooted <- midpoint(tree)
reordered_tree <- reorder(midpoint_rooted)
p <- ggtree(midpoint_rooted) + geom_tiplab(size=3) 
   # + geom_text(aes(label=label, subset=.(!isTip)))
dd <- read.table("tables/esxMatrix.txt", sep="\t", stringsAsFactor=F)
dd[dd=="no"] <- ""
p <- gheatmap(p, dd, offset = 1.0, width = 0.5, colnames_position="top", font.size=3) + 
    theme(legend.position="none") +
    scale_fill_manual(values=c("black")) + geom_treescale()
p
## can't find a way to rotate the labels, so may have to be done in illustrator
## can't find a way to add bootstrap values the way i'd like

