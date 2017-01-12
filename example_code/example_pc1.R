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

# Full ordering can be loaded here:
data("Vignette_sort")
CCLE.samp.sort <- Vignette_sort[[1]]

CCLE.pc1 <- PC1VecFun(top.gem = CCLE.mito,
                      seed.sort = CCLE.samp.sort,
                      n = 10)

CCLE.cor.vec <-  CVEval(gem.part = CCLE.mito,
                            gem.all = CCLE_small,
                            seed = CCLE.seed,
                            splits = 10)

CCLE.bic <- ThresholdBic(cor.vec = CCLE.cor.vec,sort.order = CCLE.samp.sort,
                         pc1 = as.numeric(CCLE.pc1))

CCLE.pc1 <- PC1Align(gem = CCLE_small, pc1 = CCLE.pc1,
                     cor.vec = CCLE.cor.vec ,
                     sort.order = CCLE.samp.sort,
                     bic =CCLE.bic)

CCLE.fork <- ForkClassifier(CCLE.pc1, samp.num = length(CCLE.bic[[2]]))

