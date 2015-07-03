CalcCorScoreVec <- function(ordered.list, gem, seed.size = 10, mc = FALSE,
                            num.cores = NULL){
  if(mc == TRUE){
    require(parallel)
  }
  
  order.vec <- c(seed.size:length(ordered.list))
  ol1 <- lapply(order.vec, function(x) ordered.list[seq(x)])
  
  temp.fun1 <- function(x) CorScoreCalc(gem,x)
  if(mc!=TRUE){
    cor.score.vec <- sapply(ol1, FUN = temp.fun1)
    return(cor.score.vec)
  } else {
    cor.list <- mclapply(ol1, FUN=temp.fun1, mc.cores = num.cores)
    return(unlist(cor.list))
  }

}
