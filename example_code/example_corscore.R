data(CCLE_data)
data(Mitochondrial_genes)

CCLE.data <- CCLE_data[,-c(1,2)]
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]

random.seed <- sample(seq(length = dim(CCLE.mito)[2]),10)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 100)


CorScoreCalc(CCLE.mito, random.seed)
CorScoreCalc(CCLE.mito, CCLE.seed)

CCLE.hicor.genes <- as.numeric(HclustGenesHiCor(CCLE.mito,
                                                CCLE.seed,
                                                cuts = 8))

CorScoreCalc(CCLE.mito[CCLE.hicor.genes,], CCLE.seed)
