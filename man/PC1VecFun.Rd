% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PC1VecFun.R
\name{PC1VecFun}
\alias{PC1VecFun}
\title{Calculate PC1 vector of found pattern}
\usage{
PC1VecFun(top.gem, seed.sort, n)
}
\arguments{
\item{top.gem}{Gene expression matrix containing only highly correlating genes}

\item{seed.sort}{Ordering of samples according to strength of correlation}

\item{n}{Number of samples to use in calculation of PC1}
}
\value{
PC1 value for each sample
}
\description{
The correlations found between the chosen geneset in a subset of samples
can be summarised by looking at the first principal component (PC1) 
using principal coponent analysis (PCA).
}
\details{
\code{PC1VecFun()} takes a gene expression matrix and the sample ordering
and fits a PC1 value to all the samples based on a PCA analysis done on
the first n samples.
}
\examples{
data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- (row.names(CCLE_small) \%in\% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

set.seed(102)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 100,
                      messages = 1000)

CCLE.sort <- SampleSort(gem = CCLE.mito,seed = CCLE.seed,sort.length = 11)

# Full ordering are in Vignette_sort in sysdata.rda
CCLE.samp.sort <- MCbiclust:::Vignette_sort[[1]]

CCLE.pc1 <- PC1VecFun(top.gem = CCLE.mito,
                      seed.sort = CCLE.samp.sort,
                      n = 10)

CCLE.cor.vec <-  CVEval(gem.part = CCLE.mito,
                            gem.all = CCLE_small,
                            seed = CCLE.seed,
                            splits = 10)

CCLE.bic <- ThresholdBic(cor.vec = CCLE.cor.vec,sort.order = CCLE.samp.sort,
                         pc1 = as.numeric(CCLE.pc1))

CCLE.pc1 <- PC1Align(gem = CCLE_small, pc1 = CCLE.pc1,
                     cor.vec = CCLE.cor.vec ,
                     sort.order = CCLE.samp.sort,
                     bic =CCLE.bic)

CCLE.fork <- ForkClassifier(CCLE.pc1, samp.num = length(CCLE.bic[[2]]))

}
