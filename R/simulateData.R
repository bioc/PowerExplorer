# Simulation
# INPUT
# para0.list - a list of two mean vectors of genes in two case conditions
# para1.list - a list of two theta vectors of genes in two case conditions
# simNumRep - the number of replicates of each group for simulation
# ST - simulation times
# saveResultData - whether to save the simulated data into *.RData file
# OUTPUT
# simulatedData - a list with two list objects: nullData and alterData both
# contains count data, Wald statistics, mean and dispersion of each gene
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 19th, Feb, 2017
  simulateData <- function(paraMatrix,
                           dataType=c("RNASeq", "Proteomics"),
                           isLogTransformed=FALSE,
                           enableROTS=FALSE,
                           simNumRep,
                           ST=100,
                           minLFC = 0,
                           showProcess = FALSE,
                           saveResultData=FALSE){

  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))

  numEntry <- nrow(paraMatrix)
  entryList <- rownames(paraMatrix)
  if(is.null(entryList))
    stop("Error: gene/protein names not found from input paramter matrix.")
  if(enableROTS)
    optPara <- attributes(paraMatrix)$optPara

  simulatedData <- lapply(seq_len(ST), function(simIdx){
      if(showProcess){
        percent <- simIdx*100 / ST
        text <- sprintf("%s of %s simulations", simIdx, ST)
        cat(sprintf(paste0('\r>> [%-50s] %s '),
                    paste(rep('=', percent / 2), collapse=''),
                    text))
        if (simIdx == ST)
          cat('\n')
      }
      tempMatrix <- apply(paraMatrix, 1, function(x) {
        # extract parameters: para0 - mean, para1 - dispersion/sd
        para0_0 <- x[1]; para0_1 <- x[3]
        para1_0 <- x[2]; para1_1 <- x[4]
        numMiss0 <- floor(x[5]*simNumRep[1])
        numMiss1 <- floor(x[6]*simNumRep[2])
      # simulate data for each gene/protein
      # null hypohesis: two groups should follow the same distribution
        NBdist <- dataTypeSelect & !isLogTransformed
      simData0.1 <-
        if(NBdist) {
          simNB(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                para0_0, para0_0, para1_0, para1_0)
        } else {
          simNorm(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                  para0_0, para0_0, para1_0, para1_0)
        }
      simData0.2 <-
        if(NBdist) {
          simNB(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                para0_1, para0_1, para1_1, para1_1)
        } else {
          simNorm(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                  para0_1, para0_1, para1_1, para1_1)
        }
      simData1 <-
        if(NBdist) {
          simNB(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                para0_0, para0_1, para1_0, para1_1)
        } else {
          simNorm(simNumRep[1]-numMiss0, simNumRep[2]-numMiss1,
                  para0_0, para0_1, para1_0, para1_1)
          }

      if(enableROTS){
        statistics0.1 <-
          try(ROTS.stat(simdata = simData0.1, opt.para = optPara))
        statistics0.2 <-
          try(ROTS.stat(simdata = simData0.2, opt.para = optPara))
        statistics1 <-
          try(ROTS.stat(simdata = simData1, opt.para = optPara))
      }
      else{
        statistics0.1 <-
          if(NBdist) try(WaldTest(simData0.1))
          else try(tTestGLM(simData0.1))
        statistics0.2 <-
          if(NBdist) try(WaldTest(simData0.2))
        else try(tTestGLM(simData0.2))
        statistics1 <-
          if(NBdist) try(WaldTest(simData1))
          else try(tTestGLM(simData1))
      }

      # if error occurs, test statistics remain as 0
      if (inherits(statistics0.1,"try-error")) {
        statistics0.1 <- 0
      }
      if (inherits(statistics0.2,"try-error")) {
        statistics0.2 <- 0
      }

      if (inherits(statistics1,"try-error")) {
        statistics1 <- 0
      }
      statistics0 <-
        max(statistics0.1, statistics0.2)
      # zero statistics considered as error entry
      # set both to zero to introduce penalty
      if(statistics0==0) statistics1 <- 0

      # collection all the simulated counts and statistics
      if(statistics0.1 >= statistics0.2)
        data0 <- simData0.1[,1]
      else
        data0 <- simData0.2[,1]

      data1 <- simData1[,1]
      
      res_data <- c(as.numeric(data0), rep(NA, (numMiss0+numMiss1)),
                    as.numeric(data1), rep(NA, (numMiss0+numMiss1)),
                    statistics0, statistics1)
      return(res_data)
      })
      rownames(tempMatrix) <- c(rep(c("data0", "data1"), each=sum(simNumRep)), 
                           "stat.null", "stat.alter")
      return(tempMatrix)
  })
  # record the attributes
  attr(simulatedData, "dataType") <- dataType
  attr(simulatedData, "repNum") <- simNumRep
  # Save data for later studies
  if(saveResultData) {
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    # simulated data name indicates the basic information
    filename.simulatedData <-
    sprintf("simulated_Data_numRep_%s_numSim_%s_%s.RData",
            paste(simNumRep, collapse ="_"),
            simNumRep,
            format(Sys.time(), "%H%M%S"))
    save(simulatedData, file=paste0(getwd(), "/savedRData/",
                                    filename.simulatedData))
    message(">> Data saved in directory savedRData.")
    message(paste0(">> Size: ", round(file.size(
    paste0(getwd(), "/savedRData/", 
           filename.simulatedData))/2^20, 2), " MB"))
  }
  simulatedData <<- simulatedData
  return(simulatedData)
}
