PC1VecFun <- function(top.gem,seed.sort,n){
  pca.matrix <- top.mat[,seed.sort[seq(length = n)]]
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