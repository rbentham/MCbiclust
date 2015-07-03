\name{SilhoetteHclust}
\alias{SilhoetteHclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Silhouette width hierarchical clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SilhoetteHclust(gem, seed, max.groups = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gem}{
%%     ~~Describe \code{gem} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{max.groups}{
%%     ~~Describe \code{max.groups} here~~
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
function (gem, seed, max.groups = 20) 
{
    require(cluster)
    cor.mat <- cor(t(gem[, seed]))
    cor.hclust <- hclust(dist(cor.mat))
    cor.hclust.dist <- dist(CCLE_mito_cor)
    sil.value <- seq(length = max.groups - 1)
    for (i in 2:10) {
        si2 <- silhouette(x = cutree(cor.hclust, k = i), dist = cor.hclust.dist)
        sil.value[i - 1] <- mean(si2[, 3])
    }
    plot(seq(length = max.groups - 1) + 1, sil.value, xlab = "k", 
        ylab = "Mean silhoette width")
    k1 <- which.max(sil.value) + 1
    si2 <- silhouette(x = cutree(cor.hclust, k = k1), dist = cor.hclust.dist)
    row.names(si2) <- row.names(gem)
    return(si2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
