library(treescape)

eccA <- read.nexus("trees/mrbayes/eccA.nexus.run2.t")
set.seed(123)
eccAtrees <- eccA[sample(500:length(eccA),200)]
eccB <- read.nexus("trees/mrbayes/eccB.nexus.run2.t")
eccBtrees <- eccB[sample(500:length(eccB),200)]
espI <- read.nexus("trees/mrbayes/espI.nexus.run2.t")
espItrees <- espI[sample(500:length(parA),200)]
virB4 <- read.nexus("trees/mrbayes/virB4.nexus.run2.t")
virB4trees <- virB4[sample(500:length(t4ss),200)]
tcpC <- read.nexus("trees/mrbayes/tcpC.nexus.run2.t")
tcpCtrees <- tcpC[sample(500:length(zincFinger),200)]

trees <- c(eccAtrees, eccBtrees, espItrees, virB4trees, tcpCtrees)
class(trees) <- "multiPhylo"
names(trees)[1:200] <- paste0("eccA",1:200)
names(trees)[201:400] <- paste0("eccB",1:200)
names(trees)[401:600] <- paste0("espI",1:200)
names(trees)[601:800] <- paste0("virB4",1:200)
names(trees)[801:1000] <- paste0("tcpC",1:200)

trees_rooted <- lapply(trees, function(x) root(x, resolve.root=TRUE, outgroup="ERR340530"))
class(trees_rooted) <- "multiPhylo"
scape <- treescape(trees_rooted,nf=6)

dtype <- c(rep("EccA", 200),rep("EccB",200),rep("EspI",200),rep("VirB4",200),rep("TcpC",200))

plotGroves(scape$pco, groups=dtype)
plotGroves(scape$pco, groups=dtype, yax=3)
plotGroves(scape$pco, groups=dtype, xax=2, yax=3)
