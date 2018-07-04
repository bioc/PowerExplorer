# extract distribution parameters from input data
#' @importFrom stats quantile
#' @importFrom vsn justvsn
#' @import ROTS
extParaMatrix <- function(dataMatrix, groupVec,
                          isLogTransformed=FALSE,
                          dataType=c("RNASeq","Proteomics"),
                          enableROTS=FALSE,
                          minLFC=0.5,
                          paraROTS=list(B=1000,
                                        K=NULL,
                                        paired=FALSE,
                                        a1=NULL,
                                        a2=NULL,
                                        progress=FALSE),
                          seed=123,
                          parallel=FALSE, 
                          BPPARAM=bpparam()){
  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  if(length(groupVec) != ncol(dataMatrix))
    stop("groupVec and sample number vary in length.")
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numRep <- as.numeric(table(groupVec))
  numGroup <- length(numRep)

  if(dataTypeSelect & !isLogTransformed) {
      cat("[DESeq2] Estimating distribution parameters...\n")
      # RNASeq raw counts and use DESeq2
      # construct a colData for DESeq2
      colData <- data.frame(group=groupVec, row.names=colnames(dataMatrix))
      res <- estimationByDESeq2(dataMatrix, colData, parallel=parallel, BPPARAM=BPPARAM)
      dispersion <- res$dispersion
      LFCRes <- res$LFCRes
      dataMatrix <- res$normCounts
      is.na(dataMatrix) <- !dataMatrix
  }
  else{
    # calcualte log fold changes between each two groups
    if(!dataTypeSelect & !isLogTransformed) {
      cat("[vsn] Transforming data...\n")
      dataMatrix <- vsn::justvsn(dataMatrix, verbose=FALSE)
      is.na(dataMatrix) <- !dataMatrix
    }
    LFCRes <- calLFC(dataMatrix, groupVec)
  }
  #     if(dataTypeSelect){
  #       cat("[vst] Transforming data...\n")
  #       colData <- data.frame(group=groupVec, row.names=colnames(dataMatrix))
  #       dataMatrix[is.na(dataMatrix)] <- 0
  #       dds <- DESeqDataSetFromMatrix(countData = dataMatrix,
  #                                     colData = colData,
  #                                     design = ~group)
  #       dataMatrix <- assay(DESeq2::vst(dds))
  #       is.na(dataMatrix) <- !dataMatrix
  #       LFCRes <- calLFC(dataMatrix, groupVec)
  #     }else{
  #       cat("[vsn] Transforming data...\n")
  #       dataMatrix <- vsn::justvsn(dataMatrix, verbose=FALSE)
  #       is.na(dataMatrix) <- !dataMatrix
  #     }
  #   }
  # }
  
  groups <- unique(groupVec)
  comp_index <- combn(seq_len(numGroup), 2)
  comp_index <- split(comp_index, col(comp_index))
  idx <- split(seq_along(groupVec),
               rep(seq_along(groups), table(groupVec)))
  names(comp_index) <- lapply(comp_index, function(x)
                        paste0(groups[x[1]], ".vs.",groups[x[2]]))
  paraMatrices <- lapply(comp_index, function(x){
    g1 <- x[1]; g2 <- x[2]
    comp_idx <- paste0(groups[g1], ".vs.",groups[g2])
    cat(sprintf("\n  ------ <%s> ------  \n", comp_idx))
    cat("\nLog2 Fold Change Quantiles:\n")
    show(round(quantile(abs(LFCRes[,comp_idx]),
                        probs=seq(0,1,0.1),
                        type=9),2))

    # extract data from each group
    data.case.g1 <- dataMatrix[, idx[[g1]]]
    data.case.g2 <- dataMatrix[, idx[[g2]]]

    if(enableROTS){
      if(is.null(paraROTS[["B"]]))
        paraROTS[["B"]] <- 1000
      if(is.null(paraROTS[["paired"]]))
        paraROTS[["paired"]] <- FALSE
      if(is.null(paraROTS[["progress"]]))
        paraROTS[["progress"]] <- FALSE
      if(is.null(paraROTS[["a1"]]) & is.null(paraROTS[["a1"]])){
        cat("\n[ROTS] Estimating statistics optimizing parameters...\n")
        rots.res <- suppressMessages(
          ROTS(data=dataMatrix[,c(idx[[g1]], idx[[g2]])],
               groups=rep(c(1,2), numRep[c(g1, g2)]),
               log=isLogTransformed,
               B=paraROTS[["B"]],
               K=paraROTS[["K"]],
               paired=paraROTS[["paired"]],
               a1=paraROTS[["a1"]],
               a2=paraROTS[["a2"]],
               progress=paraROTS[["progress"]],
               seed=seed))
        optPara <-  c(a1=rots.res$a1, a2=rots.res$a2)
        cat(sprintf("[ROTS] optimization parameters:\na1 = %s, a2 = %s \n",
                    optPara[1], optPara[2]))
      }
      else{
        cat("\n[ROTS skipped] optimization parameters a1 and a2 are provided.\n")
        optPara <- c(a1=paraROTS[["a1"]], a2=paraROTS[["a2"]])
      }
    }

    if(dataTypeSelect & !isLogTransformed) {
      #RNASeq-raw
      # calculate the means of all genes
      mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE)
      mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE)
      paraMatrix <- cbind(mean.g1, dispersion,
                          mean.g2, dispersion)
    }
    else{ #Proteomics / RNASeq-log
      meanSD.g1 <- apply(data.case.g1, 1, normDistMLE)
      meanSD.g2 <- apply(data.case.g2, 1, normDistMLE)
      paraMatrix <- t(rbind(meanSD.g1, meanSD.g2))
      colnames(paraMatrix) <-
        c("mean.g1", "sd.g1", "mean.g2", "sd.g2")
    }
    # record ratio of missing values
    numMiss.g1 <- apply(data.case.g1, 1,
                        function(x) round(mean(is.na(x)),2))
    numMiss.g2 <- apply(data.case.g2, 1,
                        function(x) round(mean(is.na(x)),2))
    paraMatrix <- cbind(paraMatrix, numMiss.g1, numMiss.g2)

    # LFC filter
    if(!is.na(minLFC) && minLFC > 0) {
      # filter the entrys with log fold change below minLFC
      abs.lfc <- abs(LFCRes[,comp_idx])
      filter <- names(abs.lfc[abs.lfc >= minLFC])
      # apply filter
      paraMatrix <- paraMatrix[filter, ]
      if(dataTypeSelect) {
        cat(sprintf("\n%s of %s genes are over minLFC threshold %s.\n",
                    length(filter), length(abs.lfc), minLFC))
      }
      else{
        cat(sprintf("\n%s of %s proteins are over minLFC threshold %s.\n",
                    length(filter), length(abs.lfc), minLFC))
      }
    }
    attr(paraMatrix, "Comparison") <- paste0(groups[x[1]], ".vs.",groups[x[2]])
    attr(paraMatrix, "numRep") <- numRep[c(g1, g2)]
    if(enableROTS)
      attr(paraMatrix, "optPara") <- optPara
    return(paraMatrix)
}) 
  attr(paraMatrices, "LFCRes") <- LFCRes
  return(paraMatrices)
}
