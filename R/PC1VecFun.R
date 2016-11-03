#' Calculate PC1 of found pattern
#' 
#' @param top.gem Gene expression matrix containing only highly correlating genes
#' @param seed.sort Ordering of samples according to strength of correlation
#' @param n Number of samples to use in calculation of PC1
#' @return PC1 value for each sample
#' @example example_code/example_pc1.R
#' @export


PC1VecFun <- function(top.gem,seed.sort,n){
    pca.matrix <- top.gem[,seed.sort[seq_len(n)]]
    pca.results <- prcomp(t(pca.matrix),center=TRUE)
    pca.loadings <- pca.results$rotation
  
    hi.cor.matrix <- top.gem[,seed.sort]
  
    pc1fun1 <- function(x){
        return(lsfit(as.matrix(pca.loadings),
                     as.matrix(hi.cor.matrix[,x]))$coef[2])}
  
    pc1.vec <- vapply(seq_len(length(seed.sort)), FUN = pc1fun1, FUN.VALUE = numeric(1))
  
    return(pc1.vec)
}