#' Make sure PC1 and Correlation Vector are aligned
#' 
#' @param gem Gene expression matrix
#' @param PC1 First principal component
#' @param CV Correlation vector
#' @param bic bicluster
#' @return Aligned PC1 vector
#' @export


PC1Align <- function(gem, PC1, CV, bic){
  
  max.cv.loc <- which.max(CV)
  cor.test <- cor(gem[max.cv.loc, bic[[2]]], PC1[seq(length = bic[[2]])])
  
  if(cor.test < 0){
    PC1 <- -PC1
  }
  return(PC1)
}