#' Calculate point score for samples based on gene groups
#' 
#' @param gem 
#' @param gloc1
#' @param gloc2
#' @return Numeric vector of point score for each sample

PointScoreCalc <- function(gem,gloc1,gloc2){
    gem1 <- gem - apply(gem,MARGIN = 1,median)
    samp.len <- dim(gem1)[2]
    
    rankings.vec <- rep(0,samp.len)
    
    for(i in seq(length = samp.len)){
        for(j in seq(length = length(gloc1))){
            if(gem1[gloc1[j], i] > 0){
                rankings.vec[i] <- rankings.vec[i] + 1
            } else{
                rankings.vec[i] <- rankings.vec[i] - 1
            }
        }
        for(k in seq(length = length(gloc2))){
            if(gem1[gloc2[k], i] > 0){
                rankings.vec[i] <- rankings.vec[i] - 1
            } else{
                rankings.vec[i] <- rankings.vec[i] + 1
            }
        }
        
    }
    return(rankings.vec)
}
