#' Classification of fork status
#' 
#' @param pc1 PC1 values of samples
#' @param samp.num Number of samples in the bicluster
#' @return Classification of each sample
#' @example example_code/example_pc1.R
#' @export

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
  fork.status[group1] <- "Upper"
  fork.status[group2] <- "Lower"
  
  return(fork.status)}

