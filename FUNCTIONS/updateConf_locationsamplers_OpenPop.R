#' @title NIMBLE script to change the sampler for AC locations to a block sampler.
#'
#' @description
#' \code{updateConf_locationSamplers_OpenPop} is a function to change the sampler for AC locations
#' Change to a block sampler for the AC x and y coordinates each year
#'
#' @examples
#' MCMCconf <- updateConf_DASamplers_OpenPop(MCMCconf = MCMCconf, indicatorVar = "sxy")

updateConf_locationSamplers_OpenPop <- function(MCMCconf, locationVar = "sxy", byRow = NULL,...){
   model <- MCMCconf$model
   locationNodes <- model$expandNodeNames(locationVar)
   locationVarInfo <- model$getVarInfo(locationVar)
   if(is.null(locationVarInfo)){stop(paste0("locationVar ", locationVar, " does not exist in the model."))}
   if(locationVarInfo$nDim != 3){stop(paste0("locationVar ", locationVar, " must be a 3-dimensional model variable"))}
   locationDims <- locationVarInfo$maxs
   groups <- array(locationNodes, c(locationDims))
   groups <- apply(groups, c(1,3), list)
   MCMCconf$removeSamplers(locationVar)
   for(g in groups){
      MCMCconf$addSampler(target = g[[1]],
                          type = "RW_block",
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE,
                          ...)
   }#g
   MCMCconf
}

