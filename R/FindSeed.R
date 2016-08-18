#' Find highly correlated seed of samples for gene expression matrix
#' 
#' @param gem Gene expression matrix
#' @param seed.size Size of sample seed
#' @param iterations Number of iterations
#' @param initial.seed Initial seed used, if NULL randomly chosen
#' @return Highly correlated seed
#' @export

FindSeed <- function (gem, seed.size, iterations, initial.seed = NULL){
  
  print(c("Iteration", "Cor Score"))
  sample.list <- list()
  sample.size <- dim(gem)[2]
  if (length(initial.seed) == seed.size){
  seed <- initial.seed
  } else {
  seed <- sample(seq(length=sample.size), seed.size)
    }
  gem.t <- t(gem)
  zero.rows <- which(apply(X = gem[, seed],MARGIN = 1,FUN = sd) == 0)

  if (length(zero.rows) != 0) {
    test.cor <- cor(gem.t[seed, -zero.rows])
  }
  else {
  test.cor <- cor(gem.t[seed, ])
  }
  test.cor.score <- mean(abs(test.cor))
  
  rv <- sample(seq(length=seed.size), iterations, replace = TRUE)
  # Randomly deselect one sample, replace with randomly chosen non-selected sample
  # if new seed has higher correlation score, keep new seed, if not keep.
  for (i in seq(length = iterations)) {
    seed2 <- seed[-rv[i]]
    avoid.samples <- seed2
    seed2 <- c(seed2, sample(seq(length = sample.size)[-avoid.samples],1))
        zero.rows <- which(apply(X = gem[, seed2],MARGIN = 1,FUN = sd) == 0)
        if (length(zero.rows) != 0) {
          test.cor <- cor(gem.t[seed2, -zero.rows])
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
              pre.replace <- c(1:(taken.out-1))
              post.replace <- c((taken.out+1):seed.size)
              seed <- c(seed[pre.replace],seed2[seed.size],seed[post.replace])
              }
            test.cor.score <- test.cor.score2
           
        }
        if (i%%100 == 0) {
          print(c(i, test.cor.score))
        }
    }

      return(seed)
}
