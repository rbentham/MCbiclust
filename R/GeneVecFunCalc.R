#' Function needed for gene vector calculation
#' 
#' @param gem Gene set expression matrix
#' @param seed Seed of highly correlating samples
#' @param n Number of splits

#' @return Calculation

GeneVecFunCalc <- function(gem,seed,n){
    test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = n)
    test2 <- lapply(seq_len(n),FUN = function(x) (test1 == x))
    temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed) * 
      sqrt(sum(test2[[x]],na.rm = T))
    test3 <- sapply(seq_len(n), FUN = temp.fun)
    return(test3)
}