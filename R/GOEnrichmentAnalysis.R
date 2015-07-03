GOEnrichmentAnalysis <- function(gene.names,gene.values,sig.rate,sig.type = "Mann-Whitney"){
  data(GO_term_matrix)
  data(GO_term_genes)
  if(sig_type=="Mann-Whitney"){
    pvalues<-MannWhitneyGOTerms(Genes,GeneValues)
  }
  else  if(sig_type=="ztest"){
    pvalues<-ZTestGOTerms(Genes,GeneValues)
  } else{
    stop('sig_type must be  "Mann-Whitney" or "ztest", see help for details.')}
  
  GenelistNum <- function(y){
    return(length(which(gene.names %in% y)))
  }
  GenelistLoc <- function(y){
    return((which(gene.names %in% y)))
  }
  GenelistMean <- function(y){
    return(mean(gene.values[y]))
  }
  sig.p <- which((pvalues) < sig.rate)
  ordering.p <- order(pvalues[sig.p])
  num.genes <- as.numeric(unlist(lapply(GO_term_genes[sig.p],length)))[ordering.p]
  g.in.genelist <- as.numeric(lapply(GO_term_genes[sig.p],GenelistNum))[ordering.p]
  g.av.value <- as.numeric(lapply(lapply(GO_term_genes[sig.p],GenelistLoc),GenelistMean))[ordering.p]
  p.value <- pvalues[sig.p][ordering.p]
  
  return(cbind(GO_term_matrix[sig.p[ordering.p],],num.genes,g.in.genelist,p.value,g.av.value))
}