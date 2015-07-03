\name{GOEnrichmentAnalysis}
\alias{GOEnrichmentAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GO terms enrichment analysis.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GOEnrichmentAnalysis(gene.names, gene.values, sig.rate, sig.type = "Mann-Whitney")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.names}{
%%     ~~Describe \code{gene.names} here~~
}
  \item{gene.values}{
%%     ~~Describe \code{gene.values} here~~
}
  \item{sig.rate}{
%%     ~~Describe \code{sig.rate} here~~
}
  \item{sig.type}{
%%     ~~Describe \code{sig.type} here~~
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
function (gene.names, gene.values, sig.rate, sig.type = "Mann-Whitney") 
{
    data(GO_term_matrix)
    data(GO_term_genes)
    if (sig_type == "Mann-Whitney") {
        pvalues <- MannWhitneyGOTerms(Genes, GeneValues)
    }
    else if (sig_type == "ztest") {
        pvalues <- ZTestGOTerms(Genes, GeneValues)
    }
    else {
        stop("sig_type must be  \"Mann-Whitney\" or \"ztest\", see help for details.")
    }
    GenelistNum <- function(y) {
        return(length(which(gene.names \%in\% y)))
    }
    GenelistLoc <- function(y) {
        return((which(gene.names \%in\% y)))
    }
    GenelistMean <- function(y) {
        return(mean(gene.values[y]))
    }
    sig.p <- which((pvalues) < sig.rate)
    ordering.p <- order(pvalues[sig.p])
    num.genes <- as.numeric(unlist(lapply(GO_term_genes[sig.p], 
        length)))[ordering.p]
    g.in.genelist <- as.numeric(lapply(GO_term_genes[sig.p], 
        GenelistNum))[ordering.p]
    g.av.value <- as.numeric(lapply(lapply(GO_term_genes[sig.p], 
        GenelistLoc), GenelistMean))[ordering.p]
    p.value <- pvalues[sig.p][ordering.p]
    return(cbind(GO_term_matrix[sig.p[ordering.p], ], num.genes, 
        g.in.genelist, p.value, g.av.value))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
