#' Estimate Power Under Increasing Sample Sizes
#' @description Simlilar to \code{\link{estimatePower}},
#' power estimations are performed under multiple increasing sample sizes
#' @param inputObject a numeric raw Proteomics abundance data matrix,
#' in which rows correspond to proteins and columns correspond to samples.
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; logical; set to \code{TRUE},
#' if the input data is log transformed.
#' @param dataType "RNASeq" or "Proteomics" indictes the data type of
#' the input data matrix.
#' @param minLFC LFC threshold
#' @param rangeSimNumRep a vector of sample sizes under which power
#' will be estimated
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each protein (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param enableROTS logical; if \code{TRUE}, Reproducibility-Optimized
#' Test Statistic (ROTS) will be used as the statistical model.
#' @param paraROTS a \code{list} object containing addtional parameters 
#' passed to ROTS (if enabled), see \code{\link{ROTS}}.
#' @param parallel logical; if \code{FALSE} parallelization is disabled;
#' if \code{TRUE}, parallelize calculations using 
#' \code{\link{BiocParallel}}.
#' @param BPPARAM an optional argument object passed \code{\link{bplapply}}
#' to indicate the registered cores, if \code{parallel=TRUE}.
#' @param showProcess logical; if \code{TRUE},
#' show the detailed information of
#' each simulation, used for debugging only.
#' @param saveResultData logical; if \code{TRUE}, save the simulated data
#' into RData with name pattern
#' "simulated_Data_numRep_X_numSim_XXX_XXXXX.RData".
#' @return a list of power predictions for each sample size, grouped in
#' comparisons between each two groups
#' @seealso \code{\link{estimatePower}} estimate power based on
#' actual data
#' @import utils
#' @import SummarizedExperiment
#' @importFrom Biobase exprs
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' # Example 1: a random generated Proteomics dataset (10 DE, 100 non-DE)
#' data(exampleProteomicsData)
#' dataMatrix <- exampleProteomicsData$dataMatrix
#' groupVec <- exampleProteomicsData$groupVec
#'
#' # Run estimation
#' # Note: Simulation times(ST) is specified as 5 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' predictedPower <- predictPower(dataMatrix, groupVec,
#'                                isLogTransformed=FALSE,
#'                                dataType="Proteomics",
#'                                minLFC=0,
#'                                rangeSimNumRep=c(5, 10, 15),
#'                                alpha=0.05, ST=5, seed=123)
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 13rd, March, 2018

predictPower <- function(inputObject, groupVec,
                         isLogTransformed=FALSE,
                         dataType=c("RNASeq", "Proteomics"),
                         enableROTS=FALSE,
                         paraROTS=list(B=1000,
                                       K=NULL,
                                       paired=FALSE,
                                       a1=NULL,
                                       a2=NULL,
                                       progress=FALSE),
                         minLFC = 0.5,
                         rangeSimNumRep=NA,
                         alpha=0.05, ST=100, seed=123,
                         parallel=FALSE, 
                         BPPARAM=bpparam(),
                         showProcess=FALSE,
                         saveResultData=FALSE) {
  dataMatrix <- switch (class(inputObject),
                        RangedSummarizedExperiment = assay(inputObject),
                        SummarizedExperiment = assay(inputObject),
                        ExpressionSet = Biobase::exprs(inputObject),
                        PowerExplorerStorage = {
                          if(dataType != parameters(inputObject)[["dataType"]])
                            stop("Incorrect data type.")
                          assay(inputObject)},
                        as.matrix(inputObject)
  )
  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numGroup <- length(unique(groupVec))
  timestamp()
  cat("Sample groups:\t", paste(unique(groupVec), collapse=", "),
      "\nReplicates of prediction:\t", paste(rangeSimNumRep, collapse=", "),
      "\nNum. of simulations:\t", ST,
      "\nMin. Log Fold Change:\t", minLFC,
      "\nFalse Postive Rate:\t", alpha,
      "\nTransformed:\t", isLogTransformed,
      "\nROTS enabled:\t\t\t", enableROTS,
      "\nParallel:\t\t\t", parallel,
      "\n\n")

  # remove entries with too many zero counts
  numEntries.old <- nrow(dataMatrix)
  is.na(dataMatrix) <- !dataMatrix
  # keep rows with at least two reads in each group
  exZero <-
    apply(dataMatrix, 1, function(x){
      t_groupV <- table(factor(groupVec))
      c_zero <- table(factor(groupVec)[is.na(x)])
      c_read <- t_groupV - c_zero
      return(sum(c_read <2)!=0)
    })

  dataMatrix <- dataMatrix[!exZero, ]
  cat(sprintf("%s of %s entries are filtered due to excessive zero counts\n",
              numEntries.old - nrow(dataMatrix), numEntries.old))

  # Simulation procedure
  paraMatrices <- extParaMatrix(dataMatrix=dataMatrix,
                                groupVec=groupVec,
                                isLogTransformed=isLogTransformed,
                                enableROTS=enableROTS,
                                dataType=dataType,
                                minLFC=minLFC,
                                paraROTS=paraROTS,
                                seed=seed,
                                parallel=parallel, 
                                BPPARAM=BPPARAM)

  predictedPower <- lapply(rangeSimNumRep, function(eachRepNum) {
    cat(
      sprintf(
        '\n##--Simulation with %s replicates per group--##\n', eachRepNum))

    eachRepPower <- lapply(paraMatrices, function(eachParaMatrix) {
      cat(
        sprintf(
          "\n[repNum:%s] Simulation in process, it may take a few minutes...\n",
                  eachRepNum))
      comp_idx <- attributes(eachParaMatrix)$Comparison
      # start simulation and power estimation
      cat(sprintf("\n[repNum:%s] Power Estimation between groups %s:\n",
                  eachRepNum, comp_idx))
      simData <- 
        if(!parallel){
          simulateData(eachParaMatrix,
                       dataType=dataType,
                       enableROTS=enableROTS,
                       simNumRep=c(eachRepNum, eachRepNum),
                       minLFC=minLFC,
                       ST=ST,
                       showProcess=showProcess,
                       saveResultData=saveResultData) 
        }else{
          nCores <- BPPARAM$workers
          idxCore <- factor(sort(rep(seq_len(nCores),length=ST)))
          cat(sprintf("parallelising simulations to %s workers\n", nCores))
          do.call(c, bplapply(levels(idxCore), function(i) {
            simulateData(eachParaMatrix,
                         dataType=dataType,
                         isLogTransformed=isLogTransformed,
                         enableROTS=enableROTS,
                         simNumRep=c(eachRepNum, eachRepNum),
                         minLFC=minLFC,
                         ST=sum(idxCore==i),
                         showProcess=FALSE,
                         saveResultData=saveResultData)
          }, BPPARAM=BPPARAM))
        }

      powerest <- calPwr(simData, alpha=alpha,
                         dataType=dataType,
                         saveResultData=FALSE,
                         showOverallPower=FALSE)
      temp <- apply(dataMatrix, 1, function(x) NA)
      temp[names(powerest)] <- powerest
      return(temp)
    })
    eachRepPower <- do.call(cbind, eachRepPower)
    return(eachRepPower)
  })
  names(predictedPower) <- paste0("repNum: ", rangeSimNumRep)

  switch (class(inputObject),
          SummarizedExperiment = {
            SE <- inputObject
            },
          PowerExplorerStorage = {
            SE <- inputObject
          },
          {
            SE <- SummarizedExperiment(
              assays=list(counts=dataMatrix),
              rowData=S4Vectors::DataFrame(row.names=rownames(dataMatrix)))
            colData(SE) <-
              S4Vectors::DataFrame(sampleName=factor(colnames(dataMatrix)),
                        group=factor(groupVec),
                        row.names = colnames(dataMatrix))
          }
  )

  resObject <- new("PowerExplorerStorage", SE, 
                   groupVec=groupVec,
                   LFCRes=S4Vectors::DataFrame(attributes(paraMatrices)$LFCRes,
                                               check.names = FALSE),
                   parameters = list(
                     minLFC=minLFC,
                     alpha=alpha,
                     ST=ST,
                     dataType=dataType,
                     simRepNumber=rangeSimNumRep,
                     isLogTransformed=isLogTransformed,
                     enableROTS=enableROTS
                   ),
                   predPwr=predictedPower)
  if(saveResultData){
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.result <-
      sprintf(
        ifelse(dataTypeSelect,
               "RNASeq_Results_%s.RData",
               "Proteomics_Results_%s.RData"),
        format(Sys.time(), "%H%M%S")
      )
    savedRDataDir <- paste0(getwd(), "/savedRData/")
    save(resObject, file=paste0(savedRDataDir, filename.result))
    message(">> Results saved in savedRData directory.")
    message(paste0(">> Size: ", round(file.size(
      paste0(savedRDataDir, filename.result))/2^10, 2), " KB"))
  }
  message(sprintf("use listPredPower() to view the %s-wise power results.", 
                  ifelse(dataTypeSelect,"gene", "protein")))
  return(resObject)
}
