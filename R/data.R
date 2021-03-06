#' Expression of RNA spike-ins across an RNA-seq experiment
#'
#' A dataset containing the expression of Lexogen SIRV spike-ins across an experiment.
#' Counts are estimated using Kallisto and expressed in Transcript per Million (TPM).
#' See \url{https://www.lexogen.com/sirvs/downloads} for full annotation.
#' This dataset was obtained from the fold-change mix of Lot No. 216652830
#'
#' @format A list with the following entries:
#' \describe{
#' \item{tmp}{A numeric matrix with one row per exon and one column per sample, containing the exon TPM counts}
#' \item{ID}{A data.frame with gene IDs and exon IDs}
#' \item{annot}{A factor storing whether the sample contains SIRV mix E1 or E2}
#' }
'SIRV_data'
