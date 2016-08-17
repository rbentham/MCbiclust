#' Calculate numeric vector showing correlation score decrease with sample order
#' 
#' @param ordered.list Ordered list of the samples from SampleSort
#' @param gem Gene expression matrix
#' @param seed.size Size of the seed
#' @param mc True or False for whether to use parallel evaluation
#' @param num.cores Number of cores used for parallel evaluation
#' @return Vector of correlation scores
#' @export

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
