#' Log ratio copynumber of CCLE data
#'
#' A dataset containing the log ratio copynumber for over 20000 genes and nearly 1000 samples. 
#' Availiable from the broad institute: http://www.broadinstitute.org/ccle/data/browseData
#' Filename: CCLE_copynumber_byGene_2012-04-06.txt  
#'
#' @format A data frame with 23224 rows and 976 variables:
#' \itemize{
#'   \item geneName: HGNC approved gene name
#'   \item NumChr: Chromosome number
#'   \item txStart: Start position on chromosome
#'   \item txEnd: End position on chromosome
#'   \item X1321N1_CENTRAL_NERVOUS_SYSTEM: Log ratio copynumber on sample X1321N1_CENTRAL_NERVOUS_SYSTEM
#'   \item X143B_BONE: Log ratio copynumber on sample X143B_BONE
#'   \item ...
#' }
"CCLE_copy"

#' Clinical information for CCLE data
#'
#' A dataset containing clinical information for the CCLE samples. 
#' Availiable from the broad institute: http://www.broadinstitute.org/ccle/data/browseData
#' Filename: CCLE_sample_info_file_2012-04-06.txt 
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
"CCLE_samples"

#' Expression levels of CCLE data
#'
#' A dataset containing the gene-centric RMA-normalized mRNA expression data for nearly 20000 genes and nearly 1000 samples. 
#' Availiable from the broad institute: http://www.broadinstitute.org/ccle/data/browseData
#' Filename: CCLE_Expression_Entrez_2012-04-06.gct.gz  
#'
#' @format A data frame with 18988 rows and 969 variables:
#' \itemize{
#'   \item Name: Affy probe name
#'   \item Description: HGNC approved gene name
#'   \item X1321N1_CENTRAL_NERVOUS_SYSTEM: mRNA expression on sample X1321N1_CENTRAL_NERVOUS_SYSTEM
#'   \item X143B_BONE: mRNA expressionr on sample X143B_BONE
#'   \item ...
#' }
"CCLE_copy"

#' List of known mitochondrial genes
#'
#' A dataset from MitoCarta1.0 containing the 1023 mitochondrial genes
#' Availiable from the broad institute: http://www.broadinstitute.org/scientific-community/science/programs/metabolic-disease-program/publications/mitocarta/mitocarta-in-0
#'
#' @format A Character vector of the HGNC approved gene names:
"Mitochondrial_genes"

#' List of HGNC gene names in GO terms
#'
#'
#' @format A list of GO terms with corresponding gene names:
"GO_term_genes"

#' List of entrez gene names in GO terms
#'
#'
#' @format A list of GO terms with corresponding entrez gene names:
"GO_term_EG"

#' GO terms ID and description
#'
#' A dataset containing over 18000 GO terms with their ID number
#' 
#' @format A data frame with 18826 rows and 2 variables:
#'  \itemize{
#'   \item GOID: GO ID number
#'   \item TERM: GO Term description
#' }
"GO_term_matrix"
