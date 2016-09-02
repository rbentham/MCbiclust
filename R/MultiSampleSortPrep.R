#' Preparation for SampleSort for multiple patterns
#' 
#' @param gem Gene expression matrix
#' @param av.corvec List of average correlation vector 
#' @param top.genes.num Number of the top genes in correlation vector to use for sorting samples
#' @param groups List showing what runs belong to which correlation vector group
#' @param initial.seeds List of sample seeds from all runs
#' @return Order of samples by strength to correlation pattern
#' @example example_code/example_sil.R
#' @export

MultiSampleSortPrep <- function(gem, av.corvec, top.genes.num,
                                groups, initial.seeds){
  
    t.fun <- function(x){
        return(order(abs(x),decreasing = TRUE)[seq(length = top.genes.num)])
    }
  
    top.genes <- lapply(X = av.corvec,
                        FUN = t.fun)
  
    top.seed.score <- lapply(groups,
                             FUN = function(x) seq(length = length(x)))
  
    top.seed.fun1 <- function(y) sapply(seq(length = length(groups[[y]])),
                                        FUN = function(x){
                                            l1 <- top.genes[[y]]
                                            a <- groups[[y]]
                                            l2 <- initial.seeds[a][[x]]
                                            b <- as.matrix(gem)[l1,l2]
                                            return(mean(abs(cor(t(b)))))
                                            })
  
    top.seed.score <- lapply(seq(length = length(groups)),FUN = top.seed.fun1)
  
  
    t.fun2 <- function(x){
        a <- initial.seeds[groups[[x]]][[which.max(top.seed.score[[x]])[1]]]
        return(a)
    }
    top.seed <- lapply(seq(length = length(av.corvec)),
                     FUN = t.fun2)
    top.mat <- lapply(top.genes, FUN = function(x) as.matrix(gem)[x,])
  
    return(list(top.mat, top.seed))
}