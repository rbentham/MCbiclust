#' Sort all samples by strength to correlation pattern
#' 
#' @param gem 
#' @param seed
#' @param num.cores
#' @param sort.length
#' @return Order of samples by strength to correlation pattern

SampleSort <- function(gem,seed,num.cores = NULL,sort.length = NULL){
  require(compiler)
  
  seed1 <- seed
  sample.size <- dim(gem)[2]
  if(length(sort.length) == 0){
    order.size <- sample.size 
  }
  else{
    order.size <- sort.length
  }
  
  seq.vec <- seq(length = (order.size - length(seed)))  
  tcv.max <- seq(length = length(seq.vec) + 1)
  tcv.max[1] <- CorScoreCalc(gem, seed)
  temp.fun1 <- function(x) CorScoreCalc(gem,x)
for(j in seq.vec){
  next.seed <- seq(length = sample.size)[-seed1]
  len1 <- length(next.seed)
  multi.core.list <- lapply(seq(length = len1), function(x)c(seed1, next.seed[x]))
  
  if(length(num.cores)==0){
    test.cor.values <- unlist(bplapply(multi.core.list, FUN=temp.fun1,
                                       BPPARAM = MulticoreParam()))
  }else{
    test.cor.values <- unlist(bplapply(multi.core.list, FUN=temp.fun1,
                                       BPPARAM = MulticoreParam(workers = num.cores)))}
  
  tcv.max[j + 1] <- max(test.cor.values)
  seed1 <- c(seed1, next.seed[which(test.cor.values==tcv.max[j + 1])[1]])
  if(j %% 10 == 0){
    print(c(tcv.max[j+1], seed1[length(seed1)],j))}
}  

return(seed1)
}
