sourceCpp("src/CorScoreCalcAltLoop_t.cpp")

# Test time is 221
# Default 607
# Toby version is 233

SampleSortFast <- function(gem,seed,num.cores = 1,sort.length = NULL){
  message("Sort Length \t Cor Score")
  
  gem_t <- t(gem)
  
  seed1 <- seed
  sample.size <- dim(gem)[2]
  if(length(sort.length) == 0){
    order.size <- sample.size 
  }else{
    order.size <- sort.length
  }
  
  # Construct initial array to store gene sum, gene squared sum and pairwise gene sum
  seed.unary.sum = colSums(gem_t[seed,])
  seed.sq.sum = colSums(gem_t[seed,]^2)
  seed.pair.sum = t(gem_t[seed,]) %*% gem_t[seed,]
  
  sample.sq <- list()
  sample.pair <- list()
  for (sample in seq_len(nrow(gem_t))){
    sample.sq[[sample]] <- gem_t[sample,]^2
    sample.pair[[sample]] <- gem_t[sample,] %*% t(gem_t[sample,])
  }
  
  seq.vec <- seq_len(order.size - length(seed))  
  tcv.max <- seq_len(length(seq.vec) + 1)
  tcv.max[1] <- MCbiclust:::CorScoreCalc_t(gem_t, seed)
  
  if(length(num.cores)==0){
    param = bpstart(SnowParam())
  }else{
    param = bpstart(SnowParam(workers = num.cores))
  }
  
  # Change to for loop
  j = 1
  for(j in seq_len(length(seq.vec))){
    next.seed <- seq_len(sample.size)[-seed1]
    len1 <- length(next.seed)
    test.cor.values <- CorScoreCalcEigenAltLoop_t(gem_t,next.seed,sample.sq,sample.pair,length(seed1),
                                                  seed.unary.sum,seed.sq.sum,seed.pair.sum)
    
    tcv.max[j + 1] <- max(test.cor.values)
    seed1 <- c(seed1, next.seed[which(test.cor.values==tcv.max[j + 1])[1]])
    
    seed.unary.sum <- seed.unary.sum + gem_t[tail(seed1, n=1),]
    seed.sq.sum <- seed.sq.sum + gem_t[tail(seed1, n=1),]^2
    seed.pair.sum <- seed.pair.sum + gem_t[tail(seed1, n=1),] %*% t(gem_t[tail(seed1, n=1),])
    
    if(length(seed1) %% 10 == 0){
      message(paste(length(seed1),"\t\t", format(tcv.max[j],digits = 3)))
    }
  }
  
  bpstop(param)
  # Return only seed
  return(seed1)
}