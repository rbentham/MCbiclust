data(CCLE_data)
data(Mitochondrial_genes)

CCLE.data <- CCLE_data[,-c(1,2)]
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]

set.seed(102)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 6000,
                      messages = 1000)

CCLE.sort <- SampleSort(gem = CCLE.mito,seed = CCLE.seed,sort.length = 20)

# Full ordering can be loaded here:
CCLE.samp.sort <- read.csv("~/Documents/PhD/Biclustering Paper 2016/Vignette_samplesort.txt", sep="")[,1]


CCLE.pc1 <- PC1VecFun(top.gem = CCLE.mito,
                      seed.sort = CCLE.samp.sort,
                      n = 10)

CCLE.cor.vec <-  CVEval(gem.part = CCLE.mito,
                            gem.all = CCLE.data,
                            seed = CCLE.seed,
                            splits = 10)

CCLE.bic <- ThresholdBic(cor.vec = CCLE.cor.vec,sort.order = CCLE.samp.sort,
                         pc1 = as.numeric(CCLE.pc1))

CCLE.pc1 <- PC1Align(gem = CCLE.data, PC1 = CCLE.pc1,
                     CV = CCLE.cor.vec ,
                     sort.order = CCLE.samp.sort,
                     bic =CCLE.bic)

CCLE.fork <- ForkClassifier(CCLE.pc1,samp.num = length(CCLE.bic[[2]]))

