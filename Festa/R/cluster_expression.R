#' Detect constitutive and alternative exons
#'
#' Cluster genes according to reciprocal correlations, then iteratively cut tree (bottom up)
#' until one cluster is the most expressed or tied for expression across all samples.
#' Uses \code{\link[amap:hcluster]{hcluster}} from amap package for clustering
#'
#' @param data A data.frame or numeric matrix with one row per exon and one column per sample, containing the expression values of each exon.
#' @param exonID A data.frame with the following columns. \describe{
#' \item{geneID}{Gene ID of each exon.}
#' \item{exonID}{ID of each exon.}
#' }
#' @param signDigits An integer. Number of digits rounded from expression scores for ranking calculations.
#' @param exceptions Integer. Sets the maximum number of samples in which the constitutive exon group does not need to be the most expressed one. Defaults to 10\% of the dataset.
#' @param distMethod The distance metric used by the clustering algorithm (default: correlation), see function \code{\link[amap:hcluster]{hcluster}} from package amap for more information.
#' @param link Agglomeration method used by the clustering algorithm (default: complete), see function \code{\link[amap:hcluster]{hcluster}} from package amap for more information.
#' @param nbproc Integer. Number of subprocess for parallelization (default: 1), see function \code{\link[amap:hcluster]{hcluster}} from package amap for more information
#'
#' @return A data.frame with the following columns. \describe{
#' \item{geneID}{Gene ID of each exon.}
#' \item{clusters}{Numeric ID for each exon group detected.}
#' \item{exonID}{Exon ID.}
#' \item{splicing_category}{Gene class. Can be either spliced (isoforms detected), unspliced (no isoforms detected) or  single_expressed_exon (only one exon present in the gene).}
#' \item{constitutive}{Identifies whether the exon group is constitutive (main gene) or facultative (isoform).}
#' \item{transcriptID}{Identifier for each exon group detected. Ends in '_con' for constitutive and '_fac' for facultative exons.}
#' }
#' @examples
#'  ExonAssignment = ClusterExons(data = SIRV_data$tpm, exonID = SIRV_data$ID)
#'
#' @export
ClusterExons <- function(data = NULL, exonID  = NULL,
                         exceptions = ceiling(ncol(data)*.1), signDigits = 3,
                         distMethod = "correlation", link = "complete", nbproc = 1) {
  ExonAssTable    <- list()
  exceptions      <- exceptions/ncol(data)
  data <- cbind(exonID, as.data.frame(data))
  row.names(data) <- data$exonID
  for (gID in unique(data$geneID)){
# Subset only exons within one gene ID, removing ID columns
    Evalues <- data[which(data$geneID%in%gID),-grep(pattern = "ID",x = names(data))]
# Annotate single-exon genes
    if (nrow(Evalues)<2) {
      ExonAssTable[[gID]] <- as.data.frame(matrix(row.names(Evalues), ncol = 1))
      names(ExonAssTable[[gID]]) <- "exonID"
      ExonAssTable[[gID]]$splicing_category <- "single_expressed_exon"
      ExonAssTable[[gID]]$clusters <- 0
      ExonAssTable[[gID]]$clusterranks <- 1
      ExonAssTable[[gID]] <- ExonAssTable[[gID]][c("clusters", "exonID", "clusterranks", "splicing_category")]
    } else {
# Generate main clustering for threshold detection
      tree <- amap::hcluster(Evalues, method = distMethod, link = link, nbproc = nbproc)
      for (nExonClusters in nrow(Evalues):1){
# Evaluate relative times it each exon group is ranked as first
        Evalues$clusters <- cutree(tree, k = nExonClusters)
        clusterranks <- plyr::ddply(Evalues, plyr::.(clusters), plyr::colwise(median))
# If there is no constitutive exon group, mark as unspliced
        if (nrow(clusterranks)==1){
          ExonAssTable[[gID]] <- as.data.frame(matrix(row.names(Evalues), ncol = 1))
          names(ExonAssTable[[gID]]) <- "exonID"
          ExonAssTable[[gID]]$clusters <- 0
          ExonAssTable[[gID]]$clusterranks <- 1
          ExonAssTable[[gID]]$splicing_category <- "unspliced"
          ExonAssTable[[gID]] <- ExonAssTable[[gID]][c("clusters", "exonID", "clusterranks", "splicing_category")]
          break} else {
            clusterranks <- apply(clusterranks[,-1], 2, function(x){
              rank(-round(x, digits = signDigits), ties.method = "min", na.last = T)
            })
            row.names(clusterranks) <- unique(Evalues$cluster)
            clusterranks <- apply(clusterranks, 1, function(x){sum(x==1)/ncol(clusterranks)})
# Control that there is only one exon group which consistently ranks 1st or tied allowing for exceptions
            if (sum((clusterranks)>=(1-exceptions))==1) {
              # create table with exon subcluster assignments
              matchmaker <- Evalues[,"clusters", drop=F]
              # store subcluster ranks
              clusters <- as.data.frame(clusterranks)
              #annotate subcluster ranks with subcluster IDs
              clusters$clusters <- c(1:nrow(clusters))
              # merge exon with subcluster ID and ranks
              matchmaker <- plyr::join(matchmaker, clusters, by="clusters", type="left")
              # add exon names
              row.names(matchmaker) <- row.names(Evalues)
              matchmaker$exonID <- row.names(matchmaker)
              ExonAssTable[[gID]] <- matchmaker
              ExonAssTable[[gID]]$splicing_category <- "spliced"
              ExonAssTable[[gID]] <- ExonAssTable[[gID]][c("clusters", "exonID", "clusterranks", "splicing_category")]
              break} else {next}}}
    }
  }
  ExonAssTable <- plyr::ldply(ExonAssTable, rbind)
  names(ExonAssTable)[1] <- "geneID"
# assign constitutive/specific status
  ExonAssTable$constitutive <- ifelse((ExonAssTable$clusterranks>=(1-exceptions)),"constitutive","facultative")
# merge cluster assignments into unique IDs with designation of constitutiveness
  ExonAssTable$transcriptID <- as.factor(paste0(
    ExonAssTable$geneID, "_t",
    ExonAssTable$clusters,
    ifelse(ExonAssTable$constitutive=="constitutive", "_con","_fac"),
    sep=""))
  # code ID variables as factors
  ExonAssTable$geneID <- as.factor(ExonAssTable$geneID)
  ExonAssTable$splicing_category <- as.factor(ExonAssTable$splicing_category)
  ExonAssTable$constitutive <- as.factor(ExonAssTable$constitutive)
  ExonAssTable[,-which(colnames(ExonAssTable)=="clusterranks")]
}

#' Average expression values based on unique eigenexon IDs
#'
#' @param data A data.frame or numeric matrix with one row per exon and one column per sample, containing the expression values of each exon.
#' @param spliceID A data.frame. Output from the \code{\link{ClusterExons}} function.
#' @param splicingRatios Logical. If FALSE, expression from all entries is reported on the same scale. If TRUE, expression from splicing entries is normalized by their gene's constitutive expression score, generating splicing ratios.
#' @param NAcorrection: Logical. Applicable only if splicingRatios is TRUE. If TRUE, splicing ratios higher than 1 are set to 1 and NA/NaN/infinity values to 0. This accounts for experimental error in measurements.
#' @examples
#'  ExonAssignment = ClusterExons(data = SIRV_data$tpm, exonID = SIRV_data$ID)
#'  AverageExons(SIRV_data$tpm, ExonAssignment)

AverageExons <- function(data = NULL, spliceID = NULL, splicingRatios = F, NAcorrection = F){
  spliceID = spliceID[,c("geneID", "transcriptID", "constitutive")]
  data = cbind(spliceID, data.frame(data))
  if (splicingRatios == F) {
    out <- plyr::ddply(.data = data, .variables = plyr::.(transcriptID), plyr::numcolwise(median), na.rm = T)
    out[order(out$transcriptID),]
  } else {
    ## split into constitutives and facultatives
    ConTranscripts <- data[which(data$constitutive=="constitutive"),]
    FacTranscripts <- data[which(data$constitutive!="constitutive"),]
    ## average exon values within transcripts
    ConTranscripts <- plyr::ddply(ConTranscripts, .variables = plyr::.(geneID, transcriptID), plyr::numcolwise(median), na.rm = T)
    FacTranscripts <- plyr::ddply(FacTranscripts, .variables = plyr::.(geneID, transcriptID), plyr::numcolwise(median), na.rm = T)
    FacSplicing <- apply(FacTranscripts, 1, function(Fac){
      Spl <- as.numeric(Fac[-grep(pattern = "ID", x = names(FacTranscripts))])
      Con <- ConTranscripts[which(Fac["geneID"]==ConTranscripts$geneID),-grep(pattern = "ID", x = names(ConTranscripts))]
      Spl/Con
    })
    FacSplicing <- cbind(FacTranscripts[,c("geneID","transcriptID")], plyr::ldply(FacSplicing))
    if (NAcorr == T) {
      # set NA/NaN and infinity scores to 0, set scores greater than 1 to 1
      FacSplicing[,sapply(FacSplicing, is.numeric)] <- apply(FacSplicing[,sapply(FacSplicing, is.numeric)],c(1,2),function(x){
        as.numeric(ifelse(is.na(x), 0, x))})
      FacSplicing[,sapply(FacSplicing, is.numeric)] <- apply(FacSplicing[,sapply(FacSplicing, is.numeric)],c(1,2),function(x){
        as.numeric(ifelse(x>1, 1, x))})
    }
    out <- rbind(ConTranscripts, FacSplicing)
    out[order(out$transcriptID),]
  }
}
