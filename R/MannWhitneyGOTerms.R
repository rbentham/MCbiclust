#' Runs Mann-Whitney test for all GO terms
#' 
#' @param genes Gene names
#' @param gene.values Gene values
#' @param GO_term_genes Genes in each GO term
#' @return p-value for each GO Term

MannWhitneyGOTerms <- function(genes, gene.values,GO_term_genes){
  
    mannfun1 <- function(x){
        a <- (genes %in% GO_term_genes[[x]])
        if(sum(a,na.rm = T) > 10){
          return(wilcox.test(gene.values[a],
                      gene.values)$p.value)
        }else{
          return(NA)
        }}
        
    go.pvalues <- sapply(seq_len(length(GO_term_genes)), FUN = mannfun1)
    return((go.pvalues))
}