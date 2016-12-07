#' Method for the calculation of a correlation vector
#' 
#' @param gem.part Part of gene expression matrix only containing gene set of interest
#' @param gem.all All of gene expression matrix
#' @param seed Seed of highly correlating samples
#' @param splits Number of cuts from hierarchical clustering
#' @return Correlation vector
#' @example example_code/example_pc1.R
#' @name CVEval

NULL

#' @export
#' @rdname CVEval
CVEval <- function(gem.part, gem.all, seed, splits){
      gene.vec <- GeneVecFun(gem.part, seed, splits)
      dim1 <- dim(gem.all)[1]
      temp.fun <- function(x) return(cor(as.numeric(gem.all[x, seed]), gene.vec))
      return(vapply(seq_len(dim1), temp.fun,FUN.VALUE = numeric(1)))
}


# Calcuclates gene vector used in calculation of correlation vector

GeneVecFun <- function(gem,seed,splits){
  
  test.list <- lapply(X = seq_len(splits)[-1],
                      FUN = function(x) GeneVecFunCalc(gem, seed, x))
  
  test4 <- which.max(vapply(test.list,max,
                            FUN.VALUE = numeric(1)))
  test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = (test4 + 1))
  test2 <- lapply(seq_len((test4 + 1)),
                  FUN = function(x) (test1 == x))
  
  temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed) * 
    sqrt(length(test2[[x]]))
  
  test3 <- vapply(seq_len((test4 + 1)), FUN = temp.fun, FUN.VALUE = numeric(1))
  test5 <- test2[[which.max(test3)]]
  return(colMeans(gem[test5,seed]))
}

# Function needed for gene vector calculation

GeneVecFunCalc <- function(gem,seed,n){
  test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = n)
  test2 <- lapply(seq_len(n),FUN = function(x) (test1 == x))
  temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed) * 
    sqrt(sum(test2[[x]],na.rm = TRUE))
  test3 <- vapply(seq_len(n), FUN = temp.fun, FUN.VALUE = numeric(1))
  return(test3)
}
