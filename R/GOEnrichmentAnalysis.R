#' Calculate gene set enrichment of correlation vector using Mann-Whitney test
#' 
#' @param gene.names Names of the genes
#' @param gene.values Values associated with the genes 
#' @param sig.rate Level of significance

#' @return Significant gene sets
#' @example example_code/example_GOEnrichment.R
#' @export

GOEnrichmentAnalysis <- function(gene.names,gene.values,sig.rate){
  GO.data <-GoDataLoad()
  
  GO_term_matrix <- GO.data[[1]]
  GO_term_genes <- GO.data[[2]]

  pvalues <- MannWhitneyGOTerms(gene.names,gene.values,GO_term_genes)
  adj.pvalues <- p.adjust(pvalues)

  GenelistNum <- function(y){
    return(length(which(gene.names %in% y)))
  }
  GenelistLoc <- function(y){
    return((which(gene.names %in% y)))
  }
  GenelistMean <- function(y){
    return(mean(gene.values[y]))
  }
  sig.p <- which((adj.pvalues) < sig.rate)
  ordering.p <- order(adj.pvalues[sig.p])
  num.genes <- as.numeric(unlist(lapply(GO_term_genes[sig.p],
                                        length)))[ordering.p]
  g.in.genelist <- as.numeric(lapply(GO_term_genes[sig.p],
                                     GenelistNum))[ordering.p]
  CV.av.value <- as.numeric(lapply(lapply(GO_term_genes[sig.p],GenelistLoc),
                                   GenelistMean))[ordering.p]
  p.value <- pvalues[sig.p][ordering.p]
  adj.p.value <- adj.pvalues[sig.p][ordering.p]
  
  av.gene.values <- mean(gene.values)
  
  phenofun <- function(x) return(ifelse(CV.av.value[x] > av.gene.values,1,-1))
  phenotype <- sapply(seq(length = length(sig.p)),FUN = phenofun)
  
  return(cbind(GO_term_matrix[sig.p[ordering.p],], num.genes,
               g.in.genelist, adj.p.value, CV.av.value, phenotype))


}