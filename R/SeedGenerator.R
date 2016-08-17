#' Calculate initial seeds to minimise overlaps
#' 
#' @param seed.size Size of seed
#' @param numbers Number of seeds to be generated
#' @param sample.length Number of samples
#' @param break.num Target length of intersections with previous seed
#' @param attempts Number of attempts before best result is selected
#' @return Initial seed list
#' @export

SeedGenerator <- function(seed.size,numbers,sample.length,break.num = 0, attempts = 1000){
  sample_list3 <- list()
  sample_list3[[1]] <- sample(c(1:sample.length),seed.size)
  for(i in 2:numbers){
    candidate_list <- list()
    intersect_list <- list()
    for(j in 1:attempts){
      candidate_list[[j]] <- sample(c(1:sample.length),10)
      temp.fun1 <- function(x){return(length(intersect(x,candidate_list[[j]])))}
      intersect_list[[j]] <- max(unlist(lapply(X = sample_list3[1:(i-1)], FUN = temp.fun1)))
      if(intersect_list[[j]] <= break.num){
        break
      }
    }
    sample_list3[[i]] <- candidate_list[[which.min(unlist(intersect_list))[1]]]
  }
  return(sample_list3)
}