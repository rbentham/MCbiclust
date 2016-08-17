#' Sort all samples by strength to a known correlation pattern
#' 
#' @param gem Gene expression matrix
#' @param seed Seed for highly correlating samples
#' @param group1.loc Location of genes highly correlated to each other but anti-correlated to group2
#' @param group2.loc Location of genes highly correlated to each other but anti-correlated to group1
#' @param num.cores Number of cores used in parallel evaluation
#' @param sort.length Number of samples to be sorted
#' @return Order of samples by strength to correlation pattern
#' @export

SampleSortGroups <-
  function(gem,group1.loc,group2.loc,seed,num.cores = NULL ,sort.length = NULL){

    if(length(sort.length) == 0){
      sample.size <- dim(gem)[2]}
    else{
      sample.size <- sort.length}
    
    test.cor<-cor(t(gem[,seed]))
    test.cor1 <- test.cor[group1.loc,group1.loc]
    test.cor2 <- test.cor[group2.loc,group2.loc]
    test.cor3 <- test.cor[group1.loc,group2.loc]
    
    # Normalise to group length (corresponding to area on heatmap)
    ta1 <- mean(test.cor1)/(length(group1.loc)*length(group1.loc))
    ta2 <- mean(test.cor2)/(length(group2.loc)*length(group2.loc))
    ta3 <- (-mean(test.cor3)/(length(group1.loc)*length(group2.loc)))
    
    test.val.max<-c(1:(sample.size-length(seed)))
    test.val.max[1]<-mean(c(ta1,ta2,ta3)) 
    
    SampSortMulti<-function(mat1){
      test.cor<-cor(t(gem[,mat1]))
      test.cor1<-test.cor[group1.loc,group1.loc]
      test.cor2<-test.cor[group2.loc,group2.loc]
      test.cor3<-test.cor[group1.loc,group2.loc]
      
      # Normalise to group length (corresponding to area on heatmap)
      ta1<-mean(test.cor1)/(length(group1.loc)*length(group1.loc))
      ta2<-mean(test.cor2)/(length(group2.loc)*length(group2.loc))
      ta3<-(-mean(test.cor3)/(length(group1.loc)*length(group2.loc)))
      
      return(mean(c(ta1,ta2,ta3)))
    }
    
    for(j in 1:(sample.size-length(seed))){
      next.seed <- c(1:sample.size)[-c(seed)]
      test.val2 <- c(1:length(next.seed))

      multi.core.list <- rep(list(seed),length(next.seed))
      for(i1 in 1:length(multi.core.list)){
        multi.core.list[[i1]] <- c(multi.core.list[[i1]],next.seed[i1])
      }
      if(length(num.cores)==0){
        test.val2 <- unlist(bplapply(multi.core.list,FUN=SampSortMulti,
                                     BPPARAM = MulticoreParam()))
      }else{
        test.val2 <- unlist(bplapply(multi.core.list,FUN=SampSortMulti,
                                   BPPARAM = MulticoreParam(workers = num.cores)))}
      
      test.val.max[j+1] <- max(test.val2)
      seed <- c(seed,next.seed[which(test.val2==test.val.max[j+1])[1]])
      if(j %% 10 == 0){
        print(c(test.val.max[j+1],seed[length(seed)],j))}
    }  
    return(seed)
  }
