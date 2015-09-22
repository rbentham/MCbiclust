#' Calculate correlation score
#' 
#' @param gene.expr.matrix 
#' @param sample.vec
#' @return The correlation score

CorScoreCalc <- function(gene.expr.matrix,sample.vec){
  a <- abs(cor(t(gene.expr.matrix[,sample.vec])))
  return(sum(a)/length(a))
}