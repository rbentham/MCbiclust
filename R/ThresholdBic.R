#' Methods for defining a bicluster
#'
#' A bicluster is the fundamental result found using MCbiclust. These
#' three functions are essential for the precise definition of these
#' biclusters.
#' 
#' \code{ThresholdBic()} takes as its main inputs the correlation vector 
#' (output of \code{CVEval()}), sample ordering (output of
#' \code{SampleSort()}), PC1 vector (output of \code{PC1VecFun}) and returns
#' a list of the genes and samples which belong to the bicluster according 
#' to a certain level of significance.
#' 
#' \code{PC1Align()} is a function used once the bicluster has been found to
#' ensure that the upper fork samples (those with higher PC1 values) correspond
#' to those samples that have genes with positive correlation vector values
#'up-regulated.
#' 
#' \code{ForkClassifier()} is a function used to classify which samples are in
#' the upper or lower fork.
#'  
#'
#' @param gem Gene expression matrix containing genes as rows and samples as columns.
#' @param cor.vec Correlation vector (output of \code{CVEval()}).
#' @param sort.order Order of samples (output of \code{SampleSort()}).
#' @param pc1 PC1 values for samples (output of \code{PC1VecFun}).
#' @param samp.sig Value between 0 and 1 determining number of samples in bicluster
#' @param bic bicluster (output of \code{ThresholdBic()})
#' @param samp.num Number of samples in the bicluster
#' @return Defined bicluster
#' @example example_code/example_pc1.R
#' @name ThresholdBic

NULL

#' @export
#' @rdname ThresholdBic

ThresholdBic <- function(cor.vec,sort.order,pc1,samp.sig = 0){
  
    cor.vec.kmeans <- kmeans(abs(cor.vec),centers = 2)
    genes.group1.loc <- (cor.vec.kmeans$cluster == 1)
    genes.group2.loc <- (cor.vec.kmeans$cluster == 2)
  
    if(mean(abs(cor.vec)[genes.group1.loc]) > 
       mean(abs(cor.vec)[genes.group2.loc])){
        bic.genes <- genes.group1.loc
    }else{
        bic.genes <- genes.group2.loc
    }
  
    pc1.min <- quantile(rev(pc1)[seq_len(ceiling(length(pc1) / 10))],
                        probs = 0 + (samp.sig/2))
    pc1.max <- quantile(rev(pc1)[seq_len(ceiling(length(pc1) / 10))],
                        probs = 1 - (samp.sig/2))
    first.no.samp <- which(pc1 > pc1.min & pc1 < pc1.max)[1]
  
    if(length(first.no.samp) > 0){
        bic.samps <- sort.order[seq_len(first.no.samp -1)]
    }else{
        bic.samps <- NA
    }

    return(list(bic.genes, bic.samps))
}


#' @export
#' @rdname ThresholdBic
PC1Align <- function(gem, pc1, cor.vec, sort.order, bic){
  
  max.cv.loc <- which.max(cor.vec)
  number.samps <- seq_len(length(bic[[2]]))
  gem1 <- gem[max.cv.loc, sort.order[number.samps]]
  cor.test <- cor(as.numeric(gem1), pc1[number.samps])
  
  if(cor.test < 0){
    pc1 <- -pc1
  }
  return(pc1)
}


#' @export
#' @rdname ThresholdBic
ForkClassifier <- function(pc1,samp.num){
  k.run1 <- kmeans(pc1[seq_len(samp.num)],2)
  g1 <- (k.run1[[1]] == 1)
  g2 <- (k.run1[[1]] == 2)
  if(mean(pc1[g1]) > mean(pc1[g2])){
    group1 <- g1
    group2 <- g2
  }else{
    group1 <- g2
    group2 <- g1
  }
  fork.status <- rep("None",length(pc1))
  fork.status[seq_len(samp.num)][group1] <- "Upper"
  fork.status[seq_len(samp.num)][group2] <- "Lower"
  
  return(fork.status)}

