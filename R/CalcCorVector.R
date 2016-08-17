#' Calculate correlation vector
#' 
#' @param gene.vec Single vector containing average expression levels
#' @param gem Gene expression matrix
#' @return Correlation vector
#' @export

CalcCorVector <- function(gene.vec, gem){
    dim1 <- dim(gem)[1]
    temp.fun <- function(x) return(cor(as.numeric(gem[x, ]), gene.vec))
    return(sapply(seq(length = dim1), temp.fun))
}