#' Calculate correlation score
#' 
#' The standard method to calculate the correlation score used to judge biclusters
#' in MCbiclust
#' 
#' @param gene.expr.matrix Gene expression matrix with genes as rows and samples as columns
#' @param sample.vec Vector of samples
#' @return The correlation score
#' @example example_code/example_corscore.R
#' @export

#' @export
#' @rdname CorScoreCalc

CorScoreCalc <- function(gene.expr.matrix,sample.vec){
      a <- abs(cor(t(gene.expr.matrix[,sample.vec])))
      return(sum(a,na.rm = TRUE)/length(a))
}


CorScoreCalc_t <- function(gene.expr.matrix,sample.vec){
  a <- abs(WGCNA::cor((gene.expr.matrix[sample.vec,])))
  return(sum(a,na.rm = TRUE)/length(a))
}