#' PowerExplorer Object
#' @description An extended \code{SummarizedExpriment} object to 
#' contain input dataMatrix, grouping information, estimated power, 
#' predicted power, fold change estimates and other estimation
#' parameters.
#' @export
PowerExplorerStorage <- setClass("PowerExplorerStorage",
         contains="SummarizedExperiment",
         representation=representation(
           groupVec="character",
           estPwr="DataFrame",
           predPwr="list",
           LFCRes="DataFrame",
           parameters="list"
         )
)

#' show method for PowerExplorerStorage
#' @describeIn show method for PowerExplorerStorage objects
#' @param object a PowerExplorerStorage object as input
#' @return a summary of input PowerExplorerStorage object
#' @examples
#' data(exampleObject)
#' show(exampleObject)
setMethod("show", "PowerExplorerStorage", function(object) {
  cat("##--Parameters--##\n")
  params <- parameters(object)
  cat("-dataType:", params[["dataType"]], "\n")
  cat(sprintf("-original repNum: %s\n-Comparison groups: %s\n",
              paste(as.numeric(table(groupVec(object))), collapse = ", "),
              paste(unique(groupVec(object)), collapse=", ")))
  cat(
    sprintf("-False positive rate: %s\n-LFC threshold: %s\n-Simulations: %s\n",
              params[["alpha"]], params[["minLFC"]], params[["ST"]]))
  cat("\n##--Log2 Fold Change Range--##\n")
  LFC <-
    apply(LFCRes(object), 2, function(x) round(range(abs(x), na.rm=TRUE), 2))
  rownames(LFC) <- c("minLFC", "maxLFC")
  print(LFC)
  cat("\n##--Average Estimated Power--##\n")
  if(nrow(estPwr(object)) == 0){
    cat("NA\n")
  }else{
    avgEstPwr <- apply(estPwr(object), 2, function(x) mean(x, na.rm=TRUE))
    print(avgEstPwr)
  }

  cat("\n##--Average Predicted Power--##\n")
  if(length(predPwr(object)) == 0){
    cat("NA\n")
  }else{
    avgPredPwr <- lapply(predPwr(object), function(x)
      apply(x, 2, function(y) mean(y, na.rm=TRUE)))
    print(avgPredPwr)
            }
          }
)

setGeneric("groupVec", function(object) {
  standardGeneric("groupVec")
})

setMethod("groupVec", "PowerExplorerStorage", function(object) {
  return(object@groupVec)
  })

setGeneric("groupVec<-", function(object, value) {
  standardGeneric("groupVec<-")
})

setMethod("groupVec<-", "PowerExplorerStorage", function(object, value) {
  object@groupVec <- value
  return(object)
})


setGeneric("estPwr", function(object) {
  standardGeneric("estPwr")
})

setMethod("estPwr", "PowerExplorerStorage", function(object) {
  return(object@estPwr)
})

setGeneric("estPwr<-", function(object, value) {
  standardGeneric("estPwr<-")
})

setMethod("estPwr<-", "PowerExplorerStorage", function(object, value) {
  object@estPwr <- value
  return(object)
})

setGeneric("predPwr", function(object) {
  standardGeneric("predPwr")
})

setMethod("predPwr", "PowerExplorerStorage", function(object) {
  return(object@predPwr)
})

setGeneric("predPwr<-", function(object, value) {
  standardGeneric("predPwr<-")
})

setMethod("predPwr<-", "PowerExplorerStorage", function(object, value) {
  object@predPwr <- value
  return(object)
})

setGeneric("LFCRes", function(object) {
  standardGeneric("LFCRes")
})

setMethod("LFCRes", "PowerExplorerStorage", function(object) {
  return(object@LFCRes)
})

setGeneric("LFCRes<-", function(object, value) {
  standardGeneric("LFCRes<-")
})

setMethod("LFCRes<-", "PowerExplorerStorage", function(object, value) {
  object@LFCRes <- value
  return(object)
})

setGeneric("parameters", function(object) {
  standardGeneric("parameters")
})

setMethod("parameters", "PowerExplorerStorage", function(object) {
  return(object@parameters)
})

setGeneric("parameters<-", function(object, value) {
  standardGeneric("parameters<-")
})

setMethod("parameters<-", "PowerExplorerStorage", function(object, value) {
  object@parameters <- value
  return(object)
})

