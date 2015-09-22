#' Classification of fork status
#' 
#' @param pc1
#' @param samp.num

#' @return Classification of each sample

ForkClassifier <- function(pc1,samp.num){
  k.run1 <- kmeans(pc1[seq(length = samp.num)],2)
  g1 <- which(k.run1[[1]] == 1)
  g2 <- which(k.run1[[1]] == 2)
  if(mean(pc1[g1]) > mean(pc1[g2])){
    group1 <- g1
    group2 <- g2
  }else{
    group1 <- g2
    group2 <- g1
  }
  return(list(Upper = group1, Lower = group2))}

