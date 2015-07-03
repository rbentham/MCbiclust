\name{CalcCorScoreVec}
\alias{CalcCorScoreVec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Correlation Score Vector}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalcCorScoreVec(ordered.list, gem, seed.size = 10, mc = FALSE, num.cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordered.list}{
%%     ~~Describe \code{ordered.list} here~~
}
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{seed.size}{
%%     ~~Describe \code{seed.size} here~~
}
  \item{mc}{
%%     ~~Describe \code{mc} here~~
}
  \item{num.cores}{
%%     ~~Describe \code{num.cores} here~~
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
function (ordered.list, gem, seed.size = 10, mc = FALSE, num.cores = NULL) 
{
    if (mc == TRUE) {
        require(parallel)
    }
    order.vec <- c(seed.size:length(ordered.list))
    ol1 <- lapply(order.vec, function(x) ordered.list[seq(x)])
    temp.fun1 <- function(x) CorScoreCalc(gem, x)
    if (mc != TRUE) {
        cor.score.vec <- sapply(ol1, FUN = temp.fun1)
        return(cor.score.vec)
    }
    else {
        cor.list <- mclapply(ol1, FUN = temp.fun1, mc.cores = num.cores)
        return(unlist(cor.list))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
