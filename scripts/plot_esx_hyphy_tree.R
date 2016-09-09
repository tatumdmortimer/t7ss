library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)
library(viridis)

tree <- read.tree("trees/RAxML_bipartitionsBranchLabels.esx_combine_guidance")
midpoint_rooted <- midpoint(tree)
reordered_tree <- reorder(midpoint_rooted)
hyphyData <- read.table("hyphy/guidance/esx_absrel", sep=",", header=TRUE)
hyphyResults <- read.table("hyphy/nodeLabels.txt")
hyphyResults$omega <- hyphyData$OmegaOver1
hyphyResults$p <- hyphyData$p_Holm
hyphyResults$omega[hyphyResults$p > 0.05] <- 0
hyphyResults$node <- NA
tipnode <- seq_along(reordered_tree$tip.label)
names(tipnode) <- reordered_tree$tip.label
hyphyResults$node <- tipnode[hyphyResults$V1]
i <- is.na(hyphyResults$node)
hyphyResults$node[i] = hyphyResults$V1[i]
p <- ggtree(reordered_tree, aes(color=log10(omega))) %<+% hyphyResults %>% flip(94,73) +
    geom_tiplab(size=3) +
    geom_cladelabel(node=96, label="ESX-1", align=T, offset=1) +
    geom_cladelabel(node=129, label="ESX-2", align=T, offset=1) +
    geom_cladelabel(node=106, label="ESX-3", align=T, offset=1) +
    geom_cladelabel(node=78, label="ESX-4", align=T, offset=1) +
    geom_cladelabel(node=122, label="ESX-5", align=T, offset=1) + 
    geom_treescale() +
    scale_color_viridis(option="plasma") +
    theme(legend.position="right")


