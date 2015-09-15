GeneVecFunCalc <- function(gem,seed,n){
  test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = n)
  test2 <- lapply(seq(length = n),FUN = function(x) which(test1 == x))
  temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed)*sqrt(length(test2[[x]]))
  test3 <- sapply(seq(length = n), FUN = temp.fun)
  return(test3)
}