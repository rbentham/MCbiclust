data(CCLE_data)
data(Mitochondrial_genes)

CCLE.data <- CCLE_data[,-c(1,2)]
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]

set.seed(101)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 100)

CCLE.cor.vec <- CVEval(gem.part = CCLE.mito,
                       gem.all = CCLE.data,
                       seed = CCLE.seed, splits = 10)

GEA <- GOEnrichmentAnalysis(gene.names = CCLE_data[,2],
                            gene.values = CCLE.cor.vec,
                            sig.rate = 0.05)

