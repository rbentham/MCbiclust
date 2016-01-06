#' Slihouette validation of correlation vector clusters 
#' 
#' @param cor.vec.mat
#' @param max.clusters
#' @param plots

#' @return The distinct clusters of correlation vectors


SilhouetteClustGroups <- function(cor.vec.mat, max.clusters, plots = FALSE){
  cor.vec.mat.len <- dim(cor.vec.mat)[1]
  cor.vec.mat.add <- dim(cor.vec.mat)[2] + 1
  cor.vec.mat2 <- cbind(cor.vec.mat, rnorm(cor.vec.mat.len, 0))
  
  cor.dist <- as.dist(1 - abs(cor(cor.vec.mat2)))
  cor.hclust <- hclust(cor.dist)
  
  library(cluster)
  sil.value <- seq(length = (max.clusters  - 1))
  for(i in 2:max.clusters){
    si2 <- silhouette(x = cutree(cor.hclust,k = i),dist = cor.dist)
    sil.value[i-1] <- mean(si2[,3])}
  
  if(plots == T) print(plot(seq(length = 19)+1, sil.value, xlab="Number of clusters",ylab="Mean silhoette width"))
  
  k1 <- which.max(sil.value) + 1
  si2 <- silhouette(x = cutree(cor.hclust,k = k1), dist = cor.dist)
  
  if(plots == T) print(plot(si2,col="red",main=""))
  
  cluster.groups <- lapply(seq(k1),
                           FUN=function(x) which(cutree(cor.hclust,k = k1) == x))
  
  cor.vec.minus.fun <- function(x){length(which(x==cor.vec.mat.add))}
  cor.vec.minus <- which(unlist(lapply(cluster.groups,FUN = cor.vec.minus.fun))==1)
  
  cluster.groups <- cluster.groups[-cor.vec.minus]
  return(cluster.groups)
}