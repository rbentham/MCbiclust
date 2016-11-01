data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- which(row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

CCLE.seed <- list()
CCLE.cor.vec <- list()

for(i in 1:5){
    set.seed(i)
    CCLE.seed[[i]] <- FindSeed(gem = CCLE.mito,
                               seed.size = 10,
                               iterations = 100,
                               messages = 100)}

for(i in 1:5){
  CCLE.cor.vec[[i]] <-  CVEval(gem.part = CCLE.mito,
                               gem.all = CCLE_small,
                               seed = CCLE.seed[[i]],
                               splits = 10)}

CCLE.cor.mat <- as.matrix(as.data.frame(CCLE.cor.vec))

CCLE.clust.groups <- SilhouetteClustGroups(cor.vec.mat = CCLE.cor.mat,
                                           plots = TRUE,
                                           max.clusters = 10)

av.corvec.fun <- function(x) rowMeans(CCLE.cor.mat[,x])
CCLE.average.corvec <- lapply(X = CCLE.clust.groups,
                              FUN = av.corvec.fun)

multi.sort.prep <- MultiSampleSortPrep(gem = CCLE_small,
                                       av.corvec = CCLE.average.corvec,
                                       top.genes.num = 750,
                                       groups =CCLE.clust.groups,
                                       initial.seeds = CCLE.seed)  

multi.sort <- list()        
for(i in seq_len(length(CCLE.clust.groups))){
    multi.sort[[i]] <- SampleSort(multi.sort.prep[[1]][[i]],
                                  seed = multi.sort.prep[[2]][[i]],
                                  sort.length = 11)
}

