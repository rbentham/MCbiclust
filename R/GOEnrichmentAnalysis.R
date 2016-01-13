#' Calculate gene set enrichment of correlation vector using Mann-Whitney test
#' 
#' @param gene.names
#' @param gene.values
#' @param sig.rate
#' @param output.type

#' @return Significant gene sets

GOEnrichmentAnalysis <- function(gene.names,gene.values,sig.rate,output.type = NULL){
  data(GO_term_matrix)
  data(GO_term_genes)

  
  if(length(output.type)!=0){
    if(output.type!="EM"){
      stop("output.type should be either NULL or `EM`")
    }
  }
  
  pvalues <- MannWhitneyGOTerms(gene.names,gene.values)
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
  num.genes <- as.numeric(unlist(lapply(GO_term_genes[sig.p],length)))[ordering.p]
  g.in.genelist <- as.numeric(lapply(GO_term_genes[sig.p],GenelistNum))[ordering.p]
  g.av.value <- as.numeric(lapply(lapply(GO_term_genes[sig.p],GenelistLoc),GenelistMean))[ordering.p]
  p.value <- pvalues[sig.p][ordering.p]
  adj.p.value <- adj.pvalues[sig.p][ordering.p]
  
  g.av.value[which(g.in.genelist == 0)] <- 0
  
  phenotype <- rep(-1,length(sig.p))
  for(i in seq(length = length(sig.p))){
    if(g.av.value[i] > 0) phenotype[i] <- 1
  }
  
  
  if(length(output.type) == 0){
    return(cbind(GO_term_matrix[sig.p[ordering.p],], num.genes,
                 g.in.genelist, adj.p.value, g.av.value))}
  if(output.type == "EM"){
    return(cbind(GO_term_matrix[sig.p[ordering.p],],
                 p.value, adj.p.value, phenotype))}


}