#' Methods for ordering samples
#' 
#' After finding an initial bicluster with \code{FindSeed()} the next step is
#' to extend the bicluster by ordering the remaining samples by how they 
#' preserve the correlation found.
#' 
#' \code{SampleSort()} is the basic function that achieves this, it takes the
#' gene expression matrix, seed of samples, and also has options for the number
#' of cores to run the method on and the number of samples to sort.
#' 
#' \code{MultiSampleSortPrep()} is a preparation function for \code{SampleSort()}
#' when MCbiclust has been run multiple times and returns a list of gene expression
#' matrices and seeds for each `distinct` bicluster found.
#' 
#' @param gem Gene expression matrix with genes as rows and samples as columns
#' @param seed Sample seed of highly correlating genes
#' @param num.cores Number of cores used in parallel evaluation
#' @param sort.length Number of samples to be sorted
#' @param av.corvec List of average correlation vector 
#' @param top.genes.num Number of the top genes in correlation vector to use for sorting samples
#' @param groups List showing what runs belong to which correlation vector group
#' @param initial.seeds List of sample seeds from all runs
#' @return Order of samples by strength to correlation pattern
#' @example example_code/example_pc1.R
#' @name SampleSort

NULL

#' @export
#' @rdname SampleSort

SampleSort <- function(gem,seed,num.cores = 1,sort.length = NULL){
    message("Sort Length \t Cor Score")
    
    gem_t <- t(gem)
  
    seed1 <- seed
    sample.size <- dim(gem)[2]
    if(length(sort.length) == 0){
        order.size <- sample.size 
    }
    else{
        order.size <- sort.length
    }
  
    seq.vec <- seq_len(order.size - length(seed))  
    tcv.max <- seq_len(length(seq.vec) + 1)
    tcv.max[1] <- CorScoreCalc_t(gem_t, seed)
    temp.fun1 <- function(x) CorScoreCalc_t(gem_t,x)
    
    if(length(num.cores)==0){
      param = bpstart(SnowParam())
    }else{
      param = bpstart(SnowParam(workers = num.cores))
    }
    
    for(j in seq.vec){
        next.seed <- seq_len(sample.size)[-seed1]
        len1 <- length(next.seed)
        multi.core.list <- lapply(seq_len(len1), function(x)c(seed1,
                                                                   next.seed[x]))
  
        test.cor.values <- unlist(bplapply(multi.core.list, FUN=temp.fun1,
                                           BPPARAM = param))

        tcv.max[j + 1] <- max(test.cor.values)
        seed1 <- c(seed1, next.seed[which(test.cor.values==tcv.max[j + 1])[1]])
        if(length(seed1) %% 10 == 0){
            message(paste(length(seed1),"\t\t", format(tcv.max[j+1],digits = 3)))}
        }  
    bpstop(param)
    return(seed1)
}


#' @export
#' @rdname SampleSort
MultiSampleSortPrep <- function(gem, av.corvec, top.genes.num,
                                groups, initial.seeds){
  
  t.fun <- function(x){
    return(order(abs(x),decreasing = TRUE)[seq_len(top.genes.num)])
  }
  
  top.genes <- lapply(X = av.corvec,
                      FUN = t.fun)
  
  top.seed.score <- lapply(groups,
                           FUN = function(x) seq_len(length(x)))
  
  top.seed.fun1 <- function(y) vapply(seq_len(sum(groups[[y]],na.rm = TRUE)),
                                      FUN = function(x){
                                        l1 <- top.genes[[y]]
                                        a <- groups[[y]]
                                        l2 <- initial.seeds[a][[x]]
                                        b <- as.matrix(gem)[l1,l2]
                                        return(mean(abs(cor(t(b)))))
                                      }, FUN.VALUE = numeric(1))
  
  top.seed.score <- lapply(seq_len(length(groups)),FUN = top.seed.fun1)
  
  
  t.fun2 <- function(x){
    a <- initial.seeds[groups[[x]]][[which.max(top.seed.score[[x]])[1]]]
    return(a)
  }
  top.seed <- lapply(seq_len(length(av.corvec)),
                     FUN = t.fun2)
  top.mat <- lapply(top.genes, FUN = function(x) as.matrix(gem)[x,])
  
  return(list(top.mat, top.seed))
}
