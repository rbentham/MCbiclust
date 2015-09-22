#' SampleSort for multiple patterns
#' 
#' @param gem 
#' @param av.corvec
#' @param top.genes.num
#' @param groups
#' @param initial.seeds
#' @param num.cores
#' @param sort.length
#' @return Order of samples by strength to correlation pattern

MultiSampleSort <- function(gem, av.corvec, top.genes.num, groups, initial.seeds,
                            num.cores,sort.length){
  
  top.genes <- lapply(X = av.corvec,
                      FUN = function(x) order(abs(x),decreasing = T)[seq(length = top.genes.num)])
  
  top.seed.score <- lapply(groups, FUN = function(x) seq(length = length(x)))
  
  for(j in seq(length = length(top.seed.score))){
    for(i in seq(length = length(top.seed.score[[j]]))){
      l1 <- top.genes[[j]]
      l2 <- initial.seeds[groups[[j]]][[i]]
      top.seed.score[[j]][i]<-mean(abs(cor(t(as.matrix(gem)[l1,l2]))))
    }
  }
  
  top.seed <- lapply(seq(length = length(average.corvec)),
                     FUN = function(x) initial.seeds[groups[[x]]][[which.max(top.seed.score[[x]])[1]]])
  top.mat <- lapply(top.genes, FUN = function(x) as.matrix(gem)[x,])
  
  multi.samp.sort <- lapply(seq(length = length(multi.clust.groups)),
                            FUN = function(x) SampleSort(top.mat[[x]],
                                                         seed = top.seed[[x]],num.cores = num.cores,
                                                         sort.length = sort.length))
  return(list(top.seed,top.genes,multi.samp.sort))
}