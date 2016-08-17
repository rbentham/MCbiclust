#' Calculate correlation vector
#' 
#' @param gem.part 
#' @param gem.all
#' @param seed
#' @param splits
#' @return Correlation vector

CVEval <- function(gem.part, gem.all, seed, splits){
  gene.vec <- GeneVecFun(gem.part, seed, splits)
  dim1 <- dim(gem.all)[1]
  temp.fun <- function(x) return(cor(as.numeric(gem.all[x, seed]), gene.vec))
  return(sapply(seq(length = dim1), temp.fun))
}