\name{FindSeed}
\alias{FindSeed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Seed that maximises the correlation score.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FindSeed(gem, seed.size, iterations, full.detail = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{seed.size}{
%%     ~~Describe \code{seed.size} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{full.detail}{
%%     ~~Describe \code{full.detail} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gem, seed.size, iterations, full.detail = FALSE) 
{
    sample.list <- list()
    sample.size <- dim(gem)[2]
    seed <- sample(seq(length = sample.size), seed.size)
    gem.t <- t(gem)
    zero.rows <- which(rowMeans(gem[, seed]) == 0)
    if (length(zero.rows) != 0) {
        test.cor <- cor(gem_t[seed, -zero.rows])
    }
    else {
        test.cor <- cor(gem.t[seed, ])
    }
    test.cor.score <- mean(abs(test.cor))
    if (Full_Detail == TRUE) {
        sample.list[[1]] <- c(0, test.cor.score, seed)
        sample.rank <- 1
    }
    rv <- sample(seq(length = seed.size), iterations, replace = TRUE)
    for (i in seq(length = iterations)) {
        seed2 <- seed[-rv[i]]
        avoid.samples <- seed2
        seed2 <- c(seed2, sample(seq(length = sample.size)[-avoid.samples], 
            1))
        zero.rows <- which(rowMeans(gem[, seed2]) == 0)
        if (length(zero.rows) != 0) {
            test.cor <- cor(gem.t[seed2, -zero.rows])
        }
        else {
            test.cor <- cor(gem.t[seed2, ])
        }
        test.cor.score2 <- mean(abs(test.cor))
        if (test.cor.score2 > test.cor.score) {
            taken.out <- rv[i]
            if (taken.out == 1) {
                seed <- c(seed2[seed.size], seed[-1])
            }
            else if (taken.out == seed.size) {
                seed <- c(seed[-seed.size], seed2[seed.size])
            }
            else {
                pre.replace <- c(1:(taken.out - 1))
                post.replace <- c((taken.out + 1):seed.size)
                seed <- c(seed[pre.replace], seed2[seed.size], 
                  seed[post.replace])
            }
            test.cor.score <- test.cor.score2
            if (Full_Detail == TRUE) {
                sample.rank <- sample.rank + 1
                sample.list[[sample.rank]] <- c(i, test.cor.score, 
                  seed)
            }
        }
        if (i\%\%100 == 0) {
            print(c(i, test.cor.scorel))
        }
    }
    if (Full_Detail == TRUE) {
        df_out <- as.data.frame(t(matrix(unlist(sample.list), 
            12)))
        colnames(df_out) <- c("Iteration", "Correlation Score", 
            paste("Seed", c(1:10), sep = ""))
        return(df_out)
    }
    else {
        return(main_subsamp3)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
