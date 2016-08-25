#' Calculate correlation score
#' 
#' @param gene.expr.matrix Gene expression matrix
#' @param sample.vec Vector of samples
#' @return The correlation score
#' @example example_code/example_corscore.R
#' @export

CorScoreCalc <- function(gene.expr.matrix,sample.vec){
  a <- abs(cor(t(gene.expr.matrix[,sample.vec])))
  return(sum(a)/length(a))
}