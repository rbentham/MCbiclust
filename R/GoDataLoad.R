#' Load data for Gene set enrichment analysis
#' 
#' @return List of GO terms and genes in each GO term

GoDataLoad <- function(){
    xx <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
    GO_term_list <- names(xx)
    GO_term_matrix <- (select(GO.db, GO_term_list, c("TERM","ONTOLOGY")))
  
    GO_gene_get<-function(y){
        return(as.character(unlist(mget(y,org.Hs.egSYMBOL))))}
  
    GO_term_genes<-lapply(X=xx,GO_gene_get)
  
    return(list(GO_term_matrix, GO_term_genes))
}