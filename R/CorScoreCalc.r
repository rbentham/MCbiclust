#' Calculate correlation score
#' 
#' @param gene.expr.matrix Gene expression matrix
#' @param sample.vec Vector of samples
#' @return The correlation score
#' @export

CorScoreCalc <- function(gene.expr.matrix,sample.vec){
  a <- abs(cor(t(gene.expr.matrix[,sample.vec])))
  return(sum(a)/length(a))
}