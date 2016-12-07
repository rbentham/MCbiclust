#' Calculate gene set enrichment of correlation vector using Mann-Whitney test
#' 
#' @param gene.names Names of the genes
#' @param gene.values Values associated with the genes 
#' @param sig.rate Level of significance

#' @return Significant gene sets
#' @example example_code/example_GOEnrichment.R
#' @name GOEnrichmentAnalysis

NULL

#' @export
#' @rdname GOEnrichmentAnalysis


GOEnrichmentAnalysis <- function(gene.names,gene.values,sig.rate){
    GO.data <-GoDataLoad()
  
    GO_term_matrix <- GO.data[[1]]
    GO_term_genes <- GO.data[[2]]

    pvalues <- MannWhitneyGOTerms(gene.names,gene.values,GO_term_genes)
    adj.pvalues <- p.adjust(pvalues)

    GenelistNum <- function(y){
        return(sum(gene.names %in% y, na.rm = TRUE))
    }
    GenelistLoc <- function(y){
        return(gene.names %in% y)
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
  
    phenofun <- function(x) return(ifelse(CV.av.value[x] > 
                                            av.gene.values,1,-1))
    phenotype <- vapply(seq_len(length(sig.p)),FUN = phenofun,FUN.VALUE = numeric(1))
  
    return(cbind(GO_term_matrix[sig.p[ordering.p],], num.genes,
                 g.in.genelist, adj.p.value, CV.av.value, phenotype))

}

# Runs Mann-Whitney test for all GO terms
MannWhitneyGOTerms <- function(genes, gene.values,GO_term_genes){
  
  mannfun1 <- function(x){
    a <- (genes %in% GO_term_genes[[x]])
    if(sum(a,na.rm = TRUE) > 10){
      return(wilcox.test(gene.values[a],
                         gene.values)$p.value)
    }else{
      return(NA)
    }}
  
  go.pvalues <- vapply(seq_len(length(GO_term_genes)), FUN = mannfun1, FUN.VALUE = numeric(1))
  return((go.pvalues))
}

# Load data for Gene set enrichment analysis
GoDataLoad <- function(){
  xx <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
  GO_term_list <- names(xx)
  GO_term_matrix <- AnnotationDbi::select(GO.db::GO.db, GO_term_list, c("TERM","ONTOLOGY"))
  
  GO_gene_get<-function(y){
    return(as.character(unlist(mget(y,org.Hs.egSYMBOL))))}
  
  GO_term_genes<-lapply(X=xx,GO_gene_get)
  
  return(list(GO_term_matrix, GO_term_genes))
}
