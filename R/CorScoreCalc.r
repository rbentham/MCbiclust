CorScoreCalc <- function(gene.expr.matrix,sample.vec){
  a <- abs(cor(t(gene.expr.matrix[,sample.vec])))
  return(sum(a)/length(a))
}