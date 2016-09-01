#' Slihouette validation of correlation vector clusters 
#' 
#' @param cor.vec.mat Correlation matrix of the correlation vectors (CVs)
#' @param max.clusters Maximum number of clusters to divide CVs into
#' @param plots True or False for whether to show silhouette plots
#' @param seed1 Value used to set random seed
#' @param rand.vec True or False for whether to add random correlation vector used for comparison
#' @return The distinct clusters of correlation vectors
#' @example example_code/example_sil.R
#' @export


SilhouetteClustGroups <- function(cor.vec.mat, max.clusters, plots = FALSE, seed1 = 100, rand.vec = TRUE){
  
  if(max.clusters >= dim(cor.vec.mat)[2]){
    max.clusters <- dim(cor.vec.mat)[2] - 1
  }
  
  if(rand.vec == TRUE){
    set.seed(seed1)
    cor.vec.mat.len <- dim(cor.vec.mat)[1]
    cor.vec.mat.add <- dim(cor.vec.mat)[2] + 1
    cor.vec.mat2 <- cbind(cor.vec.mat, rnorm(cor.vec.mat.len, 0))
  }else{
    cor.vec.mat2 <- cor.vec.mat
  }
  
  cor.dist <- as.dist(1 - abs(cor(cor.vec.mat2)))
  cor.hclust <- hclust(cor.dist)
  
  silfun1 <- function(x){
    si2 <- cluster::silhouette(x = cutree(cor.hclust,k = x),dist = cor.dist)
    return(mean(si2[,3]))
  }
  
  sil.value <- sapply(c(2:max.clusters), FUN = silfun1)
  
  if(plots == T) print(plot(seq(length = max.clusters-1)+1, sil.value, xlab="Number of clusters",ylab="Mean silhoette width"))
  
  k1 <- which.max(sil.value) + 1
  si2 <- cluster::silhouette(x = cutree(cor.hclust,k = k1), dist = cor.dist)
  
  if(plots == T) print(plot(si2,col="red",main=""))
  
  cluster.groups <- lapply(seq(k1),
                           FUN=function(x) which(cutree(cor.hclust,k = k1) == x))
  
  if(rand.vec == TRUE){
    cor.vec.minus.fun <- function(x){length(which(x==cor.vec.mat.add))}
    cor.vec.minus <- which(unlist(lapply(cluster.groups,FUN = cor.vec.minus.fun))==1)
  
    cluster.groups <- cluster.groups[-cor.vec.minus]
  }
  return(cluster.groups)
}