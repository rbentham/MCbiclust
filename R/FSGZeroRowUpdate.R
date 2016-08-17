#' FindSeedGroups function for updating rows with 0 sd
#' 
#' @param gem Gene expression matrix
#' @param seed Sample seed of highly correlating genes
#' @param group1.loc Location of genes highly correlated to each other but anti-correlated to group2
#' @param group2.loc Location of genes highly correlated to each other but anti-correlated to group1
#' @return Updated gene locations and correlation matrix

FSGZeroRowUpdate <- function(gem,seed,group1.loc,group2.loc){
  gem.t <- t(gem)
  zero.rows <- which(apply(X = gem[, seed],MARGIN = 1,FUN = sd) == 0)
  
  if(length(zero.rows) != 0){
    test.cor.mat <- cor(gem.t[seed, -zero.rows])
    
    group1.zr <- which(zero.rows <= length(group1.loc))
    group2.zr <- which(zero.rows > length(group1.loc))
    
    if(length(group1.zr) > 0){
      g1 <- length(group1.loc) - length(group1.zr)
      group1.loc.upd <- group1.loc[1:g1]}
    else{
      group1.loc.upd <- group1.loc}
    
    if(length(group2.zr) > 0){
      g2 <- length(group2.loc) - length(group2.zr)
      group2.loc.upd <- group2.loc[1:g2] - length(group1.zr)}
    else{
      group2.loc.upd <- group2.loc - length(group1.zr)}
  }
  else{
    test.cor.mat <- cor(t(gem[ ,seed]))
    group1.loc.upd <- group1.loc
    group2.loc.upd <- group2.loc}
  return(list(test.cor.mat,group1.loc.upd,group2.loc.upd))
}