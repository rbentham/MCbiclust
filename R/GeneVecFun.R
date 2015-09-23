#' Calcuclates gene vector used in calculation of correlation vector
#' 
#' @param gem
#' @param seed
#' @param splits

#' @return Average expression vector that matches pattern in seed

GeneVecFun <- function(gem,seed,splits){
  test.list <- list()
  for(i in 2:splits){
    test.list[[i-1]] <- GeneVecFunCalc(gem,seed,i)
  }
  test4 <- which.max(sapply(test.list,max))
  test.list[[test4]]
  test1 <- cutree(hclust(dist(cor(t(gem[,seed])))),k = (test4 + 1))
  test2 <- lapply(seq(length = (test4 + 1)),FUN = function(x) which(test1 == x))
  temp.fun <- function(x) CorScoreCalc(gem[test2[[x]],],seed) * sqrt(length(test2[[x]]))
  test3 <- sapply(seq(length = (test4 + 1)), FUN = temp.fun)
  test5 <- test2[[which.max(test3)]]
  return(colMeans(gem[test5,seed]))
}