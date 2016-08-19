#' Make sure PC1 and Correlation Vector are aligned
#' 
#' @param gem Gene expression matrix
#' @param PC1 First principal component
#' @param CV Correlation vector
#' @param sort.order Order of samples
#' @param bic bicluster
#' @return Aligned PC1 vector
#' @export


PC1Align <- function(gem, PC1, CV, sort.order, bic){
  
  max.cv.loc <- which.max(CV)
  number.samps <- seq(length = bic[[2]])
  gem1 <- gem[max.cv.loc, sort.order[number.samps]]
  cor.test <- cor(gem1, PC1[number.samps])
  
  if(cor.test < 0){
    PC1 <- -PC1
  }
  return(PC1)
}