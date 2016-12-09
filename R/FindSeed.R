#' Find highly correlated seed of samples for gene expression matrix
#' 
#' \code{FindSeed()} is the key function in MCbiclust. It takes a gene expression
#' matrix and by a stochastic method greedily searches for a seed of samples
#' that maximizes the correlation score of the chosen gene set.
#' 
#' Additional options allow for the search to start at a chosen seed, for instance
#' if a improvement to a known seed is desired. The result of \code{FindSeed()} is 
#' dependent on the number of iterations, with above 1000 usually providing a good
#' seed, and above 10000 an optimum seed.
#' 
#' @param gem Gene expression matrix with genes as rows and samples as columns
#' @param seed.size Size of sample seed
#' @param iterations Number of iterations
#' @param initial.seed Initial seed used, if NULL randomly chosen
#' @param messages frequency of progress messages
#' @return Highly correlated seed
#' @example example_code/example_corscore.R
#' @export

FindSeed <- function (gem, seed.size, iterations,
                      initial.seed = NULL, messages = 100){
  
    message("Iteration\tCorrelation Score")
    sample.list <- list()
    sample.size <- dim(gem)[2]
    if (length(initial.seed) == seed.size){
        seed <- initial.seed
    } else {
        seed <- sample(seq_len(sample.size), seed.size)
    }
    gem.t <- t(gem)
    zero.rows <- (apply(X = gem[, seed],MARGIN = 1,FUN = sd) == 0)

    if (sum(zero.rows,na.rm = TRUE) != 0) {
        test.cor <- cor(gem.t[seed, !zero.rows])
    }
    else {
        test.cor <- cor(gem.t[seed, ])
    }
    test.cor.score <- mean(abs(test.cor))
  
    rv <- sample(seq_len(seed.size), iterations, replace = TRUE)
  
    for (i in seq_len(iterations)) {
        seed2 <- seed[-rv[i]]
        avoid.samples <- seed2
        seed2 <- c(seed2, sample(seq_len(sample.size)[-avoid.samples],1))
        zero.rows <- (apply(X = gem[, seed2],MARGIN = 1,FUN = sd) == 0)
        if (sum(zero.rows,na.rm = TRUE) != 0) {
            test.cor <- cor(gem.t[seed2, !zero.rows])
        }
        else {
            test.cor <- cor(gem.t[seed2, ])
        }
        test.cor.score2 <- mean(abs(test.cor))
        if (test.cor.score2 > test.cor.score) {
            taken.out <- rv[i]
            if(taken.out == 1){
                seed <- c(seed2[seed.size], seed[-1])
            } else if(taken.out == seed.size){
                  seed <- c(seed[-seed.size], seed2[seed.size])
            } else{
                  pre.replace <- seq_len(taken.out-1)
                  post.replace <- seq_len(seed.size)[-seq_len(taken.out)]
                  seed <- c(seed[pre.replace],seed2[seed.size],
                            seed[post.replace])
                  }
            test.cor.score <- test.cor.score2
           
        }
        if (i%%messages == 0) {
            message(paste(i,"\t\t",format(test.cor.score,digits = 5)))
        }
    }

    return(seed)
}
