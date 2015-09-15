FindSeedGroups <-function(gem, seed.size, iterations, group1.loc, group2.loc,
                          initial.seed = NULL, full.detail = FALSE){
  gem <- gem[c(group1.loc,group2.loc),]
  sample.size <- dim(gem)[2]
  
  if (length(initial.seed) == seed.size){
    main.subsamp <- initial.seed}
  else {
    main.subsamp <- sample(seq(length = sample.size), seed.size)}
    
  gem.t<-t(gem)
  
  # If rows contain only zero/constant values correlation results in NA values
  zero.row.update <- FSGZeroRowUpdate(gem,main.subsamp,group1.loc,group2.loc)
  
  test.cor.mat <- zero.row.update[[1]]
  group1.loc.upd <- zero.row.update[[2]]
  group2.loc.upd <- zero.row.update[[3]]
  
  test.cor.mat1 <- test.cor.mat[group1.loc.upd, group1.loc.upd]
  test.cor.mat2 <- test.cor.mat[group2.loc.upd, group2.loc.upd]
  test.cor.mat3 <- test.cor.mat[group1.loc.upd, group2.loc.upd]
  
  # Normalise to group length (corresponding to area on heatmap)
  ta1 <- mean(test.cor.mat1) / (length(group1.loc) * length(group1.loc))
  ta2 <- mean(test.cor.mat2) / (length(group2.loc) * length(group2.loc))
  ta3 <- (-mean(test.cor.mat3) / (length(group1.loc) * length(group2.loc)))
  
  test.val <- mean(c(ta1,ta2,ta3))
  remove.sample <- sample(seq(length = seed.size), iterations, replace=TRUE)
  
  if(full.detail == TRUE){
    sample.list <- list()
    sample.list[[1]] <- c(0, test.val, mean(abs(test.cor.mat)), main.subsamp)
    sample.rank <- 1
  }
  
  for(i in 1:iterations){
    # remove one sample randomly
    test.replace <- sample(seq(length = sample.size)[-c(main.subsamp)],1)
    test.subsamp <- main.subsamp[-remove.sample[i]]
    # replace with a new sample
    test.subsamp <- c(test.subsamp, test.replace)
    
    # retest zero.rows
    zero.row.update <- FSGZeroRowUpdate(gem,test.subsamp,group1.loc,group2.loc)
    
    test.cor.mat <- zero.row.update[[1]]
    group1.loc.upd <- zero.row.update[[2]]
    group2.loc.upd <- zero.row.update[[3]]
              
    test.cor.mat1 <- test.cor.mat[group1.loc.upd,group1.loc.upd]
    test.cor.mat2 <- test.cor.mat[group2.loc.upd,group2.loc.upd]
    test.cor.mat3 <- test.cor.mat[group1.loc.upd,group2.loc.upd]
    
    # Normalise to group length (corresponding to area on heatmap)      
    ta1 <- mean(test.cor.mat1) / (length(group1.loc) * length(group1.loc))     
    ta2 <- mean(test.cor.mat2) / (length(group2.loc) * length(group2.loc))      
    ta3 <- (-mean(test.cor.mat3) / (length(group1.loc) * length(group2.loc)))
      
    test.val2 <- mean(c(ta1, ta2, ta3))
    
    if(test.val2 > test.val){
      main.subsamp <- test.subsamp
      test.val <- test.val2
      print.val <- mean(abs(test.cor.mat))
      if(full.detail == TRUE){
        sample.rank <- sample.rank + 1
        sample.list[[sample.rank]] <- c(i, test.val, print.val, main.subsamp)
      }
    }
    if(i %% 100 == 0){
      print(c(i,test.val,print.val))}
  }
  return(main.subsamp)
}

