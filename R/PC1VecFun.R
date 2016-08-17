#' Calculate PC1 of found pattern
#' 
#' @param top.gem Gene expression matrix containing only highly correlating genes
#' @param seed.sort Ordering of samples according to strength of correlation
#' @param n Number of samples to use in calculation of PC1

#' @return PC1 value for each sample
#' @export


PC1VecFun <- function(top.gem,seed.sort,n){
  pca.matrix <- top.gem[,seed.sort[seq(length = n)]]
  pca.results <- prcomp(t(pca.matrix),scores=TRUE,cor=TRUE,center=TRUE)
  pca.loadings <- pca.results$rotation
  
  hi.cor.matrix <- top.gem[,seed.sort]
  pc1.vec <- seq(length = length(seed.sort))
  
  for(i1 in seq(length = length(seed.sort))){
    pc1.vec[i1] <- lsfit(as.matrix(pca.loadings),
                         as.matrix(hi.cor.matrix[,i1]))$coef[2]
  }
  return(pc1.vec)
}