\name{ZTestGOTerms}
\alias{ZTestGOTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate p values for GO term enrichment using z test.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ZTestGOTerms(genes, gene.values, z.score = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{gene.values}{
%%     ~~Describe \code{gene.values} here~~
}
  \item{z.score}{
%%     ~~Describe \code{z.score} here~~
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
function (genes, gene.values, z.score = FALSE) 
{
    require(BSDA)
    go.term.genes.loc <- which(genes \%in\% unique(unlist(GO_term_genes)))
    genes.mean <- mean(gene.values[go.term.genes.loc])
    genes.sd <- sd(gene.values[go.term.genes.loc])
    go.terms.num <- length(GO_term_genes)
    p.values.go <- seq(length = go.terms.num)
    if (z.score == TRUE) 
        z.go <- seq(length = go.terms.num)
    for (i in seq(length = go.terms.num)) {
        genes.loc <- which(genes \%in\% GO_term_genes[[i]])
        if (length(genes.loc) > 2) {
            z.test.GO <- z.test(genes[genes.loc], y = NULL, mu = genes.mean, 
                sigma.x = genes.sd)
            p.values.go[i] <- z.test.GO$p.value
            if (z.score == TRUE) 
                z.go[i] <- z.test.GO$statistic
        }
        else {
            p.values.go[i] <- 1
            if (z.score == TRUE) 
                z.go[i] <- NA
            if (length(genes_loc) > 0) {
                mean.go[i] <- mean(gene.values[genes_loc])
            }
            else {
                mean.go[i] <- NA
            }
        }
    }
    if (z.score == TRUE) {
        return(z.go)
    }
    else {
        return(p.adjust(go.pvalues))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
