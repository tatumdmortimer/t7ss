library(treescape)
library(phangorn)
library("adegenet")
library("adegraphics")
library("ggplot2")
library(reshape2)

set.seed(123)
eccA <- read.tree("trees/mrbayes/core1009.newick")
eccAtrees <- eccA[sample(1:length(eccA),200)]
eccB <- read.tree("trees/mrbayes/core1014.newick")
eccBtrees <- eccB[sample(1:length(eccB),200)]
espI <- read.tree("trees/mrbayes/core1013.newick")
espItrees <- espI[sample(1:length(espI),200)]
virB4 <- read.tree("trees/mrbayes/core1008.newick")
virB4trees <- virB4[sample(1:length(virB4),200)]
tcpC <- read.tree("trees/mrbayes/core1007.newick")
tcpCtrees <- tcpC[sample(1:length(tcpC),200)]
nrdH <- read.tree("trees/mrbayes/core1003.newick")
nrdHtrees <- nrdH[sample(1:length(nrdH),200)]
virD <- read.tree("trees/mrbayes/core1005.newick")
virDtrees <- virD[sample(1:length(virD),200)]
eccE <- read.tree("trees/mrbayes/core1010.newick")
eccEtrees <- eccE[sample(1:length(eccE),200)]
mycP <- read.tree("trees/mrbayes/core1011.newick")
mycPtrees <- mycP[sample(1:length(mycP),200)]
eccC <- read.tree("trees/mrbayes/core1004.newick")
eccCtrees <- eccC[sample(1:length(eccC),200)]
hyp1 <- read.tree("trees/mrbayes/core1006.newick")
hyp1trees <- hyp1[sample(1:length(hyp1),200)]
eccD <- read.tree("trees/mrbayes/core1012.newick")
eccDtrees <- eccD[sample(1:length(eccD),200)]


trees <- c(eccAtrees, eccBtrees, espItrees, virB4trees, 
           tcpCtrees, nrdHtrees, virDtrees, eccEtrees, 
           mycPtrees, eccCtrees, hyp1trees, eccDtrees)
class(trees) <- "multiPhylo"
names(trees)[1:200] <- paste0("eccA",1:200)
names(trees)[201:400] <- paste0("eccB",1:200)
names(trees)[401:600] <- paste0("espI",1:200)
names(trees)[601:800] <- paste0("virB4",1:200)
names(trees)[801:1000] <- paste0("tcpC",1:200)
names(trees)[1001:1200] <- paste0("nrdH",1:200)
names(trees)[1201:1400] <- paste0("virD",1:200)
names(trees)[1401:1600] <- paste0("eccE",1:200)
names(trees)[1601:1800] <- paste0("mycP",1:200)
names(trees)[1801:2000] <- paste0("eccC",1:200)
names(trees)[2001:2200] <- paste0("hyp1",1:200)
names(trees)[2201:2400] <- paste0("eccD",1:200)


trees_rooted <- lapply(trees, function(x) midpoint(x))
class(trees_rooted) <- "multiPhylo"
scape <- treescape(trees_rooted,nf=13)

dtype <- c(rep("EccA", 200),rep("EccB",200),rep("EspI",200),rep("VirB4",200),
           rep("TcpC",200),rep("NrdH",200),rep("VirD",200),rep("EccE",200),
           rep("MycP",200),rep("EccC",200),rep("Hyp1",200),rep("EccD",200))

plotGroves(scape$pco, groups=dtype)
plotGroves(scape$pco, groups=dtype, yax=3)
plotGroves(scape$pco, groups=dtype, xax=2, yax=3)

#calculate differences between secretion system gene trees and NrdH
dist_mat <- multiDist(trees_rooted)
d <- as.matrix(dist_mat)
nrdH_dist <- as.vector(d[1001:1200,1:1000])
nrdH_dist <- append(nrdH_dist, as.vector(d[1001:1200,1201:2400]))
sec_dist <- as.vector(d[1:200,201:1000])
sec_dist <- append(sec_dist, as.vector(d[1:200, 1201:2400]))
sec_dist <- append(sec_dist, as.vector(d[201:400, 1:200]))
sec_dist <- append(sec_dist, as.vector(d[201:400, 401:1000]))
sec_dist <- append(sec_dist, as.vector(d[201:400, 1201:2400]))
sec_dist <- append(sec_dist, as.vector(d[401:600, 1:400]))
sec_dist <- append(sec_dist, as.vector(d[401:600, 601:1000]))
sec_dist <- append(sec_dist, as.vector(d[401:600, 1201:2400]))
sec_dist <- append(sec_dist, as.vector(d[601:800, 1:600]))
sec_dist <- append(sec_dist, as.vector(d[601:800, 1201:2400]))
sec_dist <- append(sec_dist, as.vector(d[601:800, 801:1000]))
sec_dist <- append(sec_dist, as.vector(d[801:1000, 1:800]))
sec_dist <- append(sec_dist, as.vector(d[801:1000, 1201:2400]))
sec_dist <- append(sec_dist, as.vector(d[1201:1400, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[1201:1400, 1401:2400]))
sec_dist <- append(sec_dist, as.vector(d[1401:1600, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[1401:1600, 1201:1400]))
sec_dist <- append(sec_dist, as.vector(d[1401:1600, 1601:2400]))
sec_dist <- append(sec_dist, as.vector(d[1601:1800, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[1601:1800, 1201:1600]))
sec_dist <- append(sec_dist, as.vector(d[1601:1800, 1801:2400]))
sec_dist <- append(sec_dist, as.vector(d[1801:2000, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[1801:2000, 1201:1800]))
sec_dist <- append(sec_dist, as.vector(d[1801:2000, 2001:2400]))
sec_dist <- append(sec_dist, as.vector(d[2001:2200, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[2001:2200, 1201:2000]))
sec_dist <- append(sec_dist, as.vector(d[2001:2200, 2201:2400]))
sec_dist <- append(sec_dist, as.vector(d[2201:2400, 1:1000]))
sec_dist <- append(sec_dist, as.vector(d[2201:2400, 1201:2200]))
df <- data.frame(sec_dist, nrdH_dist)
colnames(df) <- c(sec_name, nrd_name)
df_melt <- melt(df)
p <- ggplot(df_melt, aes(factor(variable), value)) + geom_boxplot() + theme_bw() + ylab("Tree Distance")


