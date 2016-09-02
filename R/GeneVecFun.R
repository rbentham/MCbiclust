#' Calcuclates gene vector used in calculation of correlation vector
#' 
#' @param gem Gene expression matrix
#' @param seed Seed of highly correlating samples
#' @param splits Maximum number of cuts from hierarchical clustering
#' @return Average expression vector that matches pattern in seed

GeneVecFun <- function(gem,seed,splits){
  
    test.list <- lapply(X = c(2:splits),
                        FUN = function(x) GeneVecFunCalc(gem, seed, x))
    
    test4 <- which.max(sapply(test.list,max))
    test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = (test4 + 1))
    test2 <- lapply(seq(length = (test4 + 1)),
                    FUN = function(x) which(test1 == x))
    
    temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed) * 
      sqrt(length(test2[[x]]))
    
    test3 <- sapply(seq(length = (test4 + 1)), FUN = temp.fun)
    test5 <- test2[[which.max(test3)]]
    return(colMeans(gem[test5,seed]))
}