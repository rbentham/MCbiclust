#' Find the most highly correlated genes using hierarchical clustering
#' 
#' @param gem 
#' @param seed
#' @param cuts
#' @return Numeric vector of most highly correlated genes

HclustGenesHiCor <- function(gem,seed,cuts){
  
  row.names(gem) <- seq(length = dim(gem)[1])
  
  gem.dend <- as.dendrogram(hclust(dist(cor(t(gem[,seed])))))
  gem.hclust <- hclust(dist(cor(t(gem[,seed]))))
  gem.cuts <- cutree(gem.hclust, k=cuts)
  
  hclust.genes.list <- lapply(seq(length = cuts), function(x) which(gem.cuts == x))
 
  temp.fun <- function(i) CorScoreCalc(gem[hclust.genes.list[[i]],],seed)
  hi.cor.values <- sapply(X=seq(length = cuts),temp.fun)
  
  cor.value <- CorScoreCalc(gem, seed) 
  
  return(unlist(hclust.genes.list[which(hi.cor.values > cor.value)]))
}
