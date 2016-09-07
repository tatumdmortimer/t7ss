library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

tree <- read.tree("trees/RAxML_bipartitionsBranchLabels.esx_combine_guidance")
midpoint_rooted <- midpoint(tree)
reordered_tree <- reorder(midpoint_rooted)
p <- ggtree(reordered_tree) %>% flip(94,73) +
    geom_tiplab(size=3) +
    geom_cladelabel(node=96, label="ESX-1", align=T, offset=1) +
    geom_cladelabel(node=129, label="ESX-2", align=T, offset=1) +
    geom_cladelabel(node=106, label="ESX-3", align=T, offset=1) +
    geom_cladelabel(node=78, label="ESX-4", align=T, offset=1) +
    geom_cladelabel(node=122, label="ESX-5", align=T, offset=1) + 
    geom_treescale()

p

## need to color branches based on hyphy results