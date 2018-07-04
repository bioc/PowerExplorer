# Estimate dipersion and LFC by DESeq2
# Author: Xu Qiao
# Created: 19th, Dec, 2017
# Last Modifed: 29th, Dec, 2017
#' @import DESeq2
#' 
estimationByDESeq2 <- function(dataMatrix, colData, 
                               parallel=FALSE, BPPARAM=bpparam()) {
  # initialize DESeqDataSet
  dataMatrix[is.na(dataMatrix)] <- 0
  mode(dataMatrix) <- "integer"
  dds <- DESeqDataSetFromMatrix(countData=dataMatrix,
                                colData=colData,
                                design=~ group)

  # Start DESeq dispersion estimations
  cat("\nEstimating NB parameters by DESeq2...\n")
  dds <- DESeq(dds,
               fitType="parametric",
               test="Wald",
               parallel=parallel,
               BPPARAM = BPPARAM,
               quiet=TRUE)
  # Plot dispersion estimates
  # saveGeneDispPlot(dds)
  # obtain Mean-Dispersion function form from the DESeq object
  dispersionVec <- dispersions(dds)
  names(dispersionVec) <-  rownames(dataMatrix)
  # retrieve fold changes
  LFCRes <- getLFC.rna(dds)
  output <- list(dispersion=dispersionVec, LFCRes=LFCRes)
  attr(output, "meanDispFunc") <- dispersionFunction(dds)
  return(list(dispersionVec=dispersionVec,
              LFCRes=LFCRes,
              normCounts=counts(dds, normalized=TRUE)))
}
