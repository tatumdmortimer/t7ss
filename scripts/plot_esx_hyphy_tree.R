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
hyphyResults$omega <- hyphyData$WtOmegaOver1
hyphyResults$p <- hyphyData$p_Holm
hyphyResults$omega[hyphyResults$p > 0.05] <- 0
hyphyResults$node <- NA
hyphyResults$V1 <- as.character(hyphyResults$V1)
tipnode <- seq_along(reordered_tree$tip.label)
names(tipnode) <- reordered_tree$tip.label
hyphyResults$node <- tipnode[hyphyResults$V1]
i <- is.na(hyphyResults$node)
hyphyResults$node[i] = hyphyResults$V1[i]
p <- ggtree(reordered_tree, aes(color=as.factor(prop))) %<+% hyphyResults %>% flip(94,73) +
    geom_tiplab(size=3) +
    geom_cladelabel(node=96, label="ESX-1", align=T, offset=1) +
    geom_cladelabel(node=129, label="ESX-2", align=T, offset=1) +
    geom_cladelabel(node=106, label="ESX-3", align=T, offset=1) +
    geom_cladelabel(node=78, label="ESX-4", align=T, offset=1) +
    geom_cladelabel(node=122, label="ESX-5", align=T, offset=1) + 
    scale_color_manual(values=c("gray", "blue", "green", "red", "yellow", "pink"),
                       labels=c("non-sig","<0.1", "0.1-0.2", "0.2-0.3","0.3-0.4",">0.4")) + theme(legend.position="right")


