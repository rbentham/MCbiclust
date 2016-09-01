#' Calculate bicluster threshold
#' 
#' @param cor.vec Correlation vector
#' @param sort.order Order of samples
#' @param pc1 PC1 values for samples
#' @param samp.sig Value between 0 and 1 determining number of samples in bicluster
#' @return Genes and samples in bicluster
#' @example example_code/example_pc1.R
#' @export

ThresholdBic <- function(cor.vec,sort.order,pc1,samp.sig = 0){
  
  cor.vec.kmeans <- kmeans(abs(cor.vec),centers = 2)
  genes.group1.loc <- which(cor.vec.kmeans$cluster == 1)
  genes.group2.loc <- which(cor.vec.kmeans$cluster == 2)
  
  if(mean(abs(cor.vec)[genes.group1.loc]) > mean(abs(cor.vec)[genes.group2.loc])){
    bic.genes <- genes.group1.loc
  }else{
    bic.genes <- genes.group2.loc
  }
  
  pc1.min <- quantile(rev(pc1)[seq(ceiling(length(pc1) / 10))],probs = 0 + (samp.sig/2))
  pc1.max <- quantile(rev(pc1)[seq(ceiling(length(pc1) / 10))],probs = 1 - (samp.sig/2))
  first.no.samp <- which(pc1 > pc1.min & pc1 < pc1.max)[1]
  
  if(length(first.no.samp) > 0){
    bic.samps <- sort.order[seq(length = (first.no.samp -1))]
  }else{
    bic.samps <- NA
  }

  return(list(bic.genes, bic.samps))
}