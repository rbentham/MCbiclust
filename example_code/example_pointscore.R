data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- (row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

set.seed(102)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 1000)

CCLE.sort <- SampleSort(gem = CCLE.mito,seed = CCLE.seed,sort.length = 11)

# Full ordering are in Vignette_sort in sysdata.rda
CCLE.samp.sort <- MCbiclust:::Vignette_sort[[1]]

CCLE.pc1 <- PC1VecFun(top.gem = CCLE.mito,
                      seed.sort = CCLE.samp.sort,
                      n = 10)

CCLE.hicor.genes <- as.numeric(HclustGenesHiCor(CCLE.mito,
                                                CCLE.seed,
                                                cuts = 8))

CCLE.cor.mat <- cor(t(CCLE.mito[CCLE.hicor.genes,CCLE.seed]))

gene.set1 <- labels(as.dendrogram(hclust(dist(CCLE.cor.mat)))[[1]])
gene.set2 <- labels(as.dendrogram(hclust(dist(CCLE.cor.mat)))[[2]])

gene.set1.loc <- which(row.names(CCLE.mito) %in% gene.set1)
gene.set2.loc <- which(row.names(CCLE.mito) %in% gene.set2)

ps.vec <- PointScoreCalc(CCLE.mito,gene.set1.loc,gene.set2.loc)

cor(ps.vec[CCLE.samp.sort], CCLE.pc1)
plot(ps.vec[CCLE.samp.sort])
plot(CCLE.pc1)
