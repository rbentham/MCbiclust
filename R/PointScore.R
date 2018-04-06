#' Calculate PointScore
#' 
#' Using two gene sets that are represented of a known bicluster (one gene set being 
#' up regulated while other gene set is down regulated), samples are scored based
#' on how well they match the known regulation of the bicluster.
#' 
#' The PointScore of a sample can be directly compared to the PC1 value. The PointScore
#' is typically used to identify samples related to the upper/lower fork  of a bicluster
#' without running the complete main MCbiclust pipeline on a dataset.
#' 
#' @param gene.expr.matrix Gene expression matrix with genes as rows and samples as columns
#' @param gene.loc1 Location of the rows containing the genes in gene set 1 within the gene expression matrix
#' @param gene.loc2 Location of the rows containing the genes in gene set 2 within the gene expression matrix
#' @return Vector of point scores for each sample in the gene expression matrix
#' @example example_code/example_pointscore.R
#' @export

PointScoreCalc <- function(gene.expr.matrix, gene.loc1, gene.loc2){
  if(dim(gene.expr.matrix)[2] == 1 || class(gene.expr.matrix) == 'numeric'){
    stop('Error: PointScoreCalc is not a single sample method')
  }
  if(dim(gene.expr.matrix)[2] <= 10){
    warning('PointScoreCalc is being run with very few samples. \n Method relies on the calculation of median value of each gene so results may not be accurate')
  }
  
  gem1 <- gene.expr.matrix - apply(gene.expr.matrix, MARGIN = 1, median, na.rm = TRUE)
  samp.len <- dim(gem1)[2]
  
  ps.vec <- rep(0, samp.len)
  
  for(i in seq(length = samp.len)){
    for(j in seq(length = length(gene.loc1))){
      if(is.na(gem1[gene.loc1[j], i])){
        next
      }
      if(gem1[gene.loc1[j], i] > 0){
        ps.vec[i] <- ps.vec[i] + 1
      } else{
        ps.vec[i] <- ps.vec[i] - 1
      }
    }
    for(k in seq(length = length(gene.loc2))){
      if(gem1[gene.loc2[k], i] > 0){
        ps.vec[i] <- ps.vec[i] - 1
      } else{
        ps.vec[i] <- ps.vec[i] + 1
      }
    }
  }
  return(ps.vec)
}
