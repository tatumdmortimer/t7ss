library(treescape)

eccA <- read.nexus("trees/mrbayes/eccA.nexus.run2.t")
set.seed(123)
eccAtrees <- eccA[sample(1:length(eccA),200)]
eccB <- read.nexus("trees/mrbayes/eccB.nexus.run2.t")
eccBtrees <- eccB[sample(1:length(eccB),200)]
parA <- read.nexus("trees/mrbayes/parA.nexus.run2.t")
parAtrees <- parA[sample(1:length(parA),200)]
t4ss <- read.nexus("trees/mrbayes/t4ss.nexus.run2.t")
t4sstrees <- t4ss[sample(1:length(t4ss),200)]
zincFinger <- read.nexus("trees/mrbayes/zincFinger.nexus.run2.t")
zincFingertrees <- zincFinger[sample(1:length(zincFinger),200)]

trees <- c(eccAtrees, eccBtrees, parAtrees, t4sstrees, zincFingertrees)
class(trees) <- "multiPhylo"
names(trees)[1:200] <- paste0("eccA",1:200)
names(trees)[201:400] <- paste0("eccB",1:200)
names(trees)[401:600] <- paste0("espI",1:200)
names(trees)[601:800] <- paste0("t4ss",1:200)
names(trees)[801:1000] <- paste0("zincFinger",1:200)

trees_rooted <- lapply(trees, function(x) root(x, resolve.root=TRUE, outgroup="ERR340530"))
class(trees_rooted) <- "multiPhylo"
scape <- treescape(trees_rooted,nf=6)

dtype <- c(rep("EccA", 200),rep("EccB",200),rep("EspI",200),rep("VirB4",200),rep("TcpC",200))

plotGroves(scape$pco, groups=dtype)
