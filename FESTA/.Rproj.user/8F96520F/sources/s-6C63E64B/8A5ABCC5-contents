#' Create test dataset for the ClusterExons function
#'
#' This function creates a sample test dataset for the ClusterExons function.
#' At the moment it the same number of exons per each gene.
#' Expression is sampled from a random binomial distribution.
#'
#' @param nGenes An integer, specifies the number of genes in the dataset.
#' @param nExons An integer, specifies the number of exons in each gene.
#' @param nSamples An integer, specifies the number of samples in the dataset.
#'
#' @return
#' A data.frame with nGenes*nExons rows and 2+nSamples columns.
#' The first 2 columns contain gene and exon ID.
#' Further columns contain log2 transformed expression scores.
#'
#' @examples
#' # Create a dataset with 10 genes composed by 3 exons, measured across 5 samples
#' genTestdata(10,3,5)
#'
#'  @export
genTestdata = function(nGenes, nExons, nSamples){
  geneID    = paste("gene",1:nGenes, sep = "")
  exonTable = merge(geneID,c(1:nExons))
  exonID    = paste(exonTable[ ,1],
                    exonTable[ ,2],
                    sep = "exon")
  exprData  = matrix(
    rbinom(n = nExons*nSamples, size = 1000, prob = .3),
    ncol = nSamples)
  exprData  = log2(exprData)
  data.frame(geneID = geneID,
             exonID = exonID,
             Evalue = exprData)
}
