#' Calculate numeric vector showing correlation score decrease with sample order
#' 
#' @param ordered.list 
#' @param gem
#' @param seed.size
#' @param mc
#' @param num.cores
#' @return Vector of correlation scores

CalcCorScoreVec <- function(ordered.list, gem, seed.size = 10, mc = FALSE,
                            num.cores = NULL){
  
  order.vec <- c(seed.size:length(ordered.list))
  ol1 <- lapply(order.vec, function(x) ordered.list[seq(x)])
  
  temp.fun1 <- function(x) CorScoreCalc(gem,x)
  if(mc!=TRUE){
    cor.score.vec <- sapply(ol1, FUN = temp.fun1)
    return(cor.score.vec)
  } else {
    if(length(num.cores) == 0){
      cor.list <- bplapply(ol1, FUN=temp.fun1,
                           BPPARAM = MulticoreParam())
    }else{
      cor.list <- bplapply(ol1, FUN=temp.fun1,
                         BPPARAM = MulticoreParam(workers = num.cores))}
    return(unlist(cor.list))
  }

}
