#' Clinical information for CCLE data
#'
#' A dataset containing clinical information for the CCLE samples. 
#'
#' @format A data frame with 967 rows and 14 variables:
#' \itemize{
#'   \item CCLE.name: Sample name identifier.
#'   \item Cell.line.primary.name: Cell line name.
#'   \item Cell.line.aliases: Any known aliases of cell line.
#'   \item Gender: Gender of patient cell line derived from.
#'   \item Site.Primary: Primary site cell line derived from.
#'   \item Histology: Histology of tumour cell line derived from.
#'   \item Hist.Subtype1: Histology subtype of tumour cell line derived from.
#'   \item Notes: Additional notes.
#'   \item Source: Source of the cell line.
#'   \item Expression.arrays: Expression array used.
#'   \item SNP.arrays: SNP array used.
#'   \item Oncomap: Oncomap mutation array used.
#'   \item Hybrid.Capture.Sequencing: Hybrid capture sequencing used.
#'   \item Name: Sample name identifier
#' }
#' @source \url{http://www.broadinstitute.org/ccle/data/browseData} Filename: CCLE_sample_info_file_2012-04-06.txt 
#' @return NA
"CCLE_samples"


#' Subset of expression levels of CCLE data
#'
#' A dataset containing the gene-centric RMA-normalized mRNA expression data for nearly 1000 genes and 500 samples taken as a random subset of the complete CCLE data. 
#' 1000 genes were selected randomly such that 500 were mitochondrial and 500 non-mitochondrial.
#'
#' @format A data frame with 1000 rows and 500 variables:
#' \itemize{
#'   \item MKN74_STOMACH: mRNA expression on sample MKN74_STOMACH
#'   \item OC316_OVARY: mRNA expressionr on sample OC316_OVARY
#'   \item ...
#' }
#'  @source \url{http://www.broadinstitute.org/ccle/data/browseData} Filename: CCLE_Expression_Entrez_2012-04-06.gct.gz
#' @return NA
"CCLE_small"


#' List of known mitochondrial genes
#'
#' A dataset from MitoCarta1.0 containing the 1023 mitochondrial genes
#' Availiable from the broad institute: http://www.broadinstitute.org/scientific-community/science/programs/metabolic-disease-program/publications/mitocarta/mitocarta-in-0
#'
#' @format A Character vector of the HGNC approved gene names
#' @source \url{https://www.broadinstitute.org/publications/broad807s} 
#' @return NA
"Mitochondrial_genes"

