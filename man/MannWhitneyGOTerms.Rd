\name{MannWhitneyGOTerms}
\alias{MannWhitneyGOTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate GO terms enrichment p value using Mann-Whitney test.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MannWhitneyGOTerms(genes, gene.values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{gene.values}{
%%     ~~Describe \code{gene.values} here~~
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
function (genes, gene.values) 
{
    go.pvalues <- seq(length = length(GO_term_genes))
    for (i in seq(length = length(go.pvalues))) {
        a <- which(genes \%in\% GO_term_genes[[i]])
        if (length(a) > 10) {
            go.pvalues[i] <- wilcox.test(gene.values[a], gene.values)$p.value
        }
        else {
            go.pvalues[i] <- NA
        }
    }
    return(p.adjust(go.pvalues))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
