data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- which(row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

CCLE.seed <- list()
CCLE.cor.vec <- list()

for(i in 1:3){
    set.seed(i)
    CCLE.seed[[i]] <- FindSeed(gem = CCLE.mito,
                               seed.size = 10,
                               iterations = 100,
                               messages = 100)}

for(i in 1:3){
    CCLE.cor.vec[[i]] <-  CVEval(gem.part = CCLE.mito,
                                 gem.all = CCLE_small,
                                 seed = CCLE.seed[[i]],
                                splits = 10)}



CCLE.cor.df <- (as.data.frame(CCLE.cor.vec))

CVPlot(cv.df = CCLE.cor.df, geneset.loc = mito.loc,
       geneset.name = "Mitochondrial",alpha1 = 0.5)
