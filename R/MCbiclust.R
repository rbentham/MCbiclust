#' MCbiclust: Massively Correlated biclustering
#'
#' MCbiclust is a R package for running massively correlating biclustering
#' analysis.MCbiclust aims to find large scale biclusters with selected 
#' features being highly correlated with each other over a subset of samples.
#'
#' The package was originally designed in order to solve a problem in
#' bioinformatics: to find biclusters representing different modes of regulation
#' of mitochondria gene expression in disease states such as breast cancer. 
#' The same methods however, can be used on any gene expression data set to
#' find biclusters of interest.
#'
#' To learn more about MCbiclust, start with the vignette:
#' \code{browseVignettes(package = "MCbiclust")}
#' @docType package
#' @name MCbiclust
#' @importFrom stats hclust kmeans sd p.adjust cutree dist as.dendrogram wilcox.test prcomp lsfit rnorm as.dist quantile
#' @importFrom BiocParallel bplapply MulticoreParam bpstart bpstop SnowParam
#' @importFrom graphics plot
#' @importFrom utils combn
#' @importFrom GO.db GO.db
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#' @importFrom AnnotationDbi select mget
#' @importFrom GGally ggpairs putPlot
#' @importFrom ggplot2 ggplot theme geom_density geom_point xlab ylab aes_
#' @importFrom scales hue_pal
#' @importFrom cluster silhouette
#' @importMethodsFrom AnnotationDbi as.list

NULL