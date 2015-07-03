\name{SeedGenerator}
\alias{SeedGenerator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate random seeds which are maximally different from each other.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SeedGenerator(seed.size, numbers, sample.length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed.size}{
%%     ~~Describe \code{seed.size} here~~
}
  \item{numbers}{
%%     ~~Describe \code{numbers} here~~
}
  \item{sample.length}{
%%     ~~Describe \code{sample.length} here~~
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
function (seed.size, numbers, sample.length) 
{
    sample_list3 <- list()
    sample_list3[[1]] <- sample(c(1:sample.length), seed.size)
    for (i in 2:numbers) {
        candidate_list <- list()
        intersect_list <- list()
        for (j in 1:100) {
            candidate_list[[j]] <- sample(c(1:sample.length), 
                10)
            temp.fun1 <- function(x) {
                return(length(intersect(x, candidate_list[[j]])))
            }
            intersect_list[[j]] <- max(unlist(lapply(X = sample_list3[1:(i - 
                1)], FUN = temp.fun1)))
            if (intersect_list[[j]] <= 1) {
                break
            }
        }
        print(j)
        sample_list3[[i]] <- candidate_list[[which.min(unlist(intersect_list))[1]]]
    }
    return(sample_list3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
