#' Calculate correlation vector
#' 
#' @param gene.vec 
#' @param gem
#' @return Correlation vector

CalcCorVector <- function(gene.vec, gem){
    dim1 <- dim(gem)[1]
    temp.fun <- function(x) return(cor(as.numeric(gem[x, ]), gene.vec))
    return(sapply(seq(length = dim1), temp.fun))
}