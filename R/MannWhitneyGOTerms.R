#' Runs Mann-Whitney test for all GO terms
#' 
#' @param genes 
#' @param gene.values
#' @return p-value for each GO Term

MannWhitneyGOTerms <- function(genes, gene.values){
  go.pvalues <- seq(length = length(GO_term_genes)) 
  for(i in seq(length = length(go.pvalues))){
    a <- which(genes %in% GO_term_genes[[i]])
    if(length(a) > 10){
      go.pvalues[i] <- wilcox.test(gene.values[a],gene.values)$p.value}
    else{
      go.pvalues[i] <- 1}
  }
  return((go.pvalues))
}