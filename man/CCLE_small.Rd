% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CCLE_small}
\alias{CCLE_small}
\title{Subset of expression levels of CCLE data}
\format{A data frame with 1000 rows and 500 variables:
\itemize{
  \item MKN74_STOMACH: mRNA expression on sample MKN74_STOMACH
  \item OC316_OVARY: mRNA expressionr on sample OC316_OVARY
  \item ...
}
 @source \url{http://www.broadinstitute.org/ccle/data/browseData} Filename: CCLE_Expression_Entrez_2012-04-06.gct.gz}
\usage{
CCLE_small
}
\value{
NA
}
\description{
A dataset containing the gene-centric RMA-normalized mRNA expression data for nearly 1000 genes and 500 samples taken as a random subset of the complete CCLE data. 
1000 genes were selected randomly such that 500 were mitochondrial and 500 non-mitochondrial.
}
\keyword{datasets}
