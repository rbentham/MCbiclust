\name{SampleSort}
\alias{SampleSort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sort samples by strength of correlation.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SampleSort(gem, seed, num.cores, sort.length = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{num.cores}{
%%     ~~Describe \code{num.cores} here~~
}
  \item{sort.length}{
%%     ~~Describe \code{sort.length} here~~
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
function (gem, seed, num.cores, sort.length = NULL) 
{
    require(parallel)
    seed1 <- seed
    sample.size <- dim(gem)[2]
    if (length(sort.length) == 0) {
        order.size <- sample.size
    }
    else {
        order.size <- sort.length
    }
    seq.vec <- seq(length = (sample.size - length(seed)))
    tcv.max <- seq(length = length(seq.vec) + 1)
    tcv.max[1] <- CorScoreCalc(gem, seed)
    for (j in seq.vec) {
        next.seed <- seq(length = sample.size)[-seed]
        len1 <- length(next.seed)
        multi.core.list <- lapply(seq(length = len1), function(x) c(seed, 
            next.seed[x]))
        temp.fun1 <- function(x) CorScoreCalc(gem, x)
        test.cor.values <- unlist(mclapply(multi.core.list, FUN = temp.fun1, 
            mc.cores = num.cores))
        tcv.max[j + 1] <- max(test.cor.values)
        seed1 <- c(seed1, next.seed[which(test.cor.values == 
            tcv.max[j + 1])[1]])
        if (j\%\%10 == 0) {
            print(c(tcv.max[j + 1], seed1[length(seed1)], j))
        }
    }
    return(seed1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
