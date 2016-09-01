#' Calculate correlation vector
#' 
#' @param gem.part Part of gene expression matrix only containing gene set of interest
#' @param gem.all All of gene expression matrix
#' @param seed Seed of highly correlating samples
#' @param splits Number of cuts from hierarchical clustering
#' @return Correlation vector
#' @example example_code/example_pc1.R
#' @export

CVEval <- function(gem.part, gem.all, seed, splits){
  gene.vec <- GeneVecFun(gem.part, seed, splits)
  dim1 <- dim(gem.all)[1]
  temp.fun <- function(x) return(cor(as.numeric(gem.all[x, seed]), gene.vec))
  return(sapply(seq(length = dim1), temp.fun))
}