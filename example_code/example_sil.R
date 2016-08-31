data(CCLE_data)
data(Mitochondrial_genes)

CCLE.data <- CCLE_data[,-c(1,2)]
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]

CCLE.seed <- list()
CCLE.cor.vec <- list()

for(i in 1:10){
set.seed(i)
CCLE.seed[[i]] <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 100)}

for(i in 1:10){
  CCLE.cor.vec[[i]] <-  CVEval(gem.part = CCLE.mito,
                                         gem.all = CCLE.data,
                                         seed = CCLE.seed[[i]],
                                         splits = 10)}

CCLE.cor.mat <- as.matrix(as.data.frame(CCLE.cor.vec))

CCLE.clust.groups <- SilhouetteClustGroups(cor.vec.mat = CCLE.cor.mat,
                      plots = T,
                      max.clusters = 10)

av.corvec.fun <- function(x) rowMeans(CCLE.cor.mat[,x])
CCLE.average.corvec <- lapply(X = CCLE.clust.groups,
                         FUN = av.corvec.fun)

multi.sort.prep <- MultiSampleSortPrep(gem = CCLE.data,av.corvec = CCLE.average.corvec,
                                       top.genes.num = 1000,
                                       groups =CCLE.clust.groups,
                                       initial.seeds =CCLE.seed)  

multi.sort <- list()        
for(i in seq(length = length(CCLE.clust.groups))){
multi.sort[[i]] <- SampleSort(multi.sort.prep[[1]][[i]],
                              seed = multi.sort.prep[[2]][[i]],
                              sort.length = 20)
}

