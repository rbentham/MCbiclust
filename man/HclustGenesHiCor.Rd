\name{HclustGenesHiCor}
\alias{HclustGenesHiCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical clustering identification of highly correlating genes.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HclustGenesHiCor(gem, seed, cuts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{cuts}{
%%     ~~Describe \code{cuts} here~~
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
function (gem, seed, cuts) 
{
    row.names(gem) <- seq(length = dim(gem)[1])
    gem.dend <- as.dendrogram(hclust(dist(cor(t(gem[, seed])))))
    gem.hclust <- hclust(dist(cor(t(gem[, seed]))))
    gem.cuts <- cutree(gem.hclust, k = cuts)
    hclust.genes.list <- lapply(seq(length = cuts), function(x) which(gem.cuts == 
        x))
    temp.fun <- function(i) CorScoreCalc(gem[hclust_samp_list[[i]], 
        ], seed)
    hi.cor.values <- sapply(X = c(1:8), temp.fun)
    cor.value <- CorScoreCalc(gem, seed)
    return(unlist(hclust.samp.list[which(hi.cor.values > cor.value)]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
