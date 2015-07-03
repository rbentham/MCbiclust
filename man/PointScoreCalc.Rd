\name{PointScoreCalc}
\alias{PointScoreCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Point Score Calculation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PointScoreCalc(gem, gloc1, gloc2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{gloc1}{
%%     ~~Describe \code{gloc1} here~~
}
  \item{gloc2}{
%%     ~~Describe \code{gloc2} here~~
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
function (gem, gloc1, gloc2) 
{
    gem1 <- gem - rowMeans(gem)
    samp.len <- dim(gem1)[2]
    rankings.vec <- rep(0, samp.len)
    for (i in seq(length = samp.len)) {
        for (j in seq(length = length(gloc1))) {
            if (gem1[gloc1[j], i] > 0) {
                rankings.vec[i] <- rankings.vec[i] + 1
            }
            else {
                rankings.vec[i] <- rankings.vec[i] - 1
            }
        }
        for (k in seq(length = length(gloc2))) {
            if (gem1[gloc2[k], i] > 0) {
                rankings.vec[i] <- rankings.vec[i] - 1
            }
            else {
                rankings.vec[i] <- rankings.vec[i] + 1
            }
        }
    }
    return(rankings.vec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
