#' Find the most highly correlated genes using hierarchical clustering
#' 
#' @param gem Gene expression matrix
#' @param seed Seed of highly correlating samples
#' @param cuts Number of groups to cut dendogram into
#' @return Numeric vector of most highly correlated genes
#' @example example_code/example_corscore.R
#' @export

HclustGenesHiCor <- function(gem,seed,cuts){
  
    row.names(gem) <- seq_len(dim(gem)[1])
  
    gem.dend <- as.dendrogram(hclust(dist(cor(t(gem[,seed])))))
    gem.hclust <- hclust(dist(cor(t(gem[,seed]))))
    gem.cuts <- cutree(gem.hclust, k=cuts)
  
    hclust.genes.list <- lapply(seq_len(cuts),
                                function(x) which(gem.cuts == x))
 
    temp.fun <- function(i) CorScoreCalc(gem[hclust.genes.list[[i]],],seed)
    hi.cor.values <- sapply(X=seq_len(cuts),temp.fun)
  
    cor.value <- CorScoreCalc(gem, seed) 
  
    return(unlist(hclust.genes.list[which(hi.cor.values > cor.value)]))
}
