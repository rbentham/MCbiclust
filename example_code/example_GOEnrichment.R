data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- which(row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

set.seed(101)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 100)

CCLE.cor.vec <- CVEval(gem.part = CCLE.mito,
                       gem.all = CCLE_small,
                       seed = CCLE.seed, splits = 10)

GEA <- GOEnrichmentAnalysis(gene.names = row.names(CCLE_small),
                            gene.values = CCLE.cor.vec,
                            sig.rate = 0.05)

