#' Find the most highly correlated genes using hierarchical clustering
#' 
#' Upon finding an initial bicluster with \code{FindSeed()} not all the genes
#' in the chosen geneset will be highly correlated to the bicluster. 
#' \code{HclustGenesHiCor()} uses the output of \code{FindSeed()} and 
#' hierarchical clustering to only select the genes that are most highly correlated
#' to the bicluster. This is achieved by cutting the dendogram produced
#' from the clustering into a set number of groups and then only selecting the
#' groups that are most highly correlated to the bicluster
#' 
#' @param gem Gene expression matrix with genes as rows and samples as columns
#' @param seed Seed of highly correlating samples
#' @param cuts Number of groups to cut dendogram into
#' @return Numeric vector of most highly correlated genes
#' @example example_code/example_corscore.R
#' @export

HclustGenesHiCor <- function(gem,seed,cuts){
  
  row.vec <- seq_len(dim(gem)[1])
  row.names(gem) <- seq_len(dim(gem)[1])
  sd.loc <- which(as.numeric(apply(gem, MARGIN = 1, sd)) == 0)
  
  if(length(sd.loc > 0)){
    gem <- gem[-sd.loc,]
    row.vec <- row.vec[-sd.loc]
  }
  
  gem.cor <- cor(t(gem[, seed]),
                 use = 'pairwise.complete.obs')
  
  gem.dend <- as.dendrogram(hclust(dist(gem.cor)))
  gem.hclust <- hclust(dist(gem.cor))
  gem.cuts <- cutree(gem.hclust, k=cuts)
  
  hclust.genes.list <- lapply(seq_len(cuts),
                              function(x) which(gem.cuts == x))
  
  temp.fun <- function(i) {
    CorScoreCalc(data.frame(gem)[hclust.genes.list[[i]],],seed)}
  
  hi.cor.values <- vapply(X=seq_len(cuts),temp.fun,
                          FUN.VALUE = numeric(1))
  
  cor.value <- CorScoreCalc(gem, seed) 
  hi.loc <- which(hi.cor.values > cor.value)
  row.vec[as.numeric(unlist(hclust.genes.list[hi.loc]))]
  
  return(row.vec[as.numeric(unlist(hclust.genes.list[hi.loc]))])
}

