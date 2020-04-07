#' @title NIMBLE script to create an Indicator NIMBLE sampler.
#'
#' @description
#' \code{sampler_nonzero_indicator} creates a sampler that contains another sampler and only uses it when an indicator is 1.
#' 
#' \code{updateConf_indicatorSamplers} modifies the current NIMBLE configuration to include the indicator sampler. 
#' The related variables will be sampled only if the inidcator is different from 0.
#' 
#' @examples
#' MCMCconf <- updateConf_indicatorSamplers(MCMCconf = MCMCconf, indicatorVar = "z", relatedVarList  = list("sxy", "p0", "sigma"))

#### 1.Indicator Sampler ####
sampler_nonzero_indicator <- nimbleFunction( contains = sampler_BASE
                                             , setup = function(model, mvSaved, target, control){
                                                regularSamplerInput <- control$sampler
                                                if(inherits(regularSamplerInput, "samplerConf"))
                                                   regularSamplerInput <- regularSamplerInput$buildSampler(model=model, mvSaved=mvSaved)
                                                if(!is.nf(regularSamplerInput))
                                                   stop(paste0("Problem building nonzero_indicator sampler for target",paste0(target, collapse = ",")))
                                                regularSamplers <- nimbleFunctionList(sampler_BASE)
                                                regularSamplers[[1]] <- regularSamplerInput
                                                indicatorNode <- control$indicator
                                             }
                                             , run = function(){if(model[[indicatorNode]] == 1) regularSamplers[[1]]$run()}
                                             , methods = list(reset = function() {regularSamplers[[1]]$reset()}))

#### 2.Update Function ####
updateConf_indicatorSamplers <- function( MCMCconf
                                          , indicatorVar
                                          , relatedVarList
                                          , byRow = NULL){
   model <- MCMCconf$model
   indicatorNodes <- model$expandNodeNames(indicatorVar)
   for(relatedVar in relatedVarList){
      relatedNodes <- model$expandNodeNames(relatedVar)
      relatedVarInfo <- model$getVarInfo(relatedVar)
      if(is.null(relatedVarInfo))
         stop(paste0("relatedVar ", relatedVar, " does not exist in the model."))
      relatedNDim <- relatedVarInfo$nDim
      if(relatedVarInfo$nDim > 2)
         stop(paste0("relatedVar ", relatedVar, " must be a 1- or 2- dimensional model variable"))
      relatedDims <- relatedVarInfo$maxs
      if(relatedNDim == 2) {
         if(is.null(byRow)) {
            byRow <- relatedDims[1] >= relatedDims[2]
         }
         groups <- matrix(relatedNodes, nrow = relatedDims[1], byrow = !byRow)
         groups <- apply(groups, 1, list)
      }
      if(relatedNDim == 1) {
         groups <- lapply(relatedNodes, list)
      }
      for(i in seq_along(groups)) {
         thisIndicatorNode <- indicatorNodes[i]
         thisRelatedNodes <- groups[[i]][[1]]
         existingSamplers <- MCMCconf$getSamplers(thisRelatedNodes)
         MCMCconf$removeSamplers(thisRelatedNodes)
         for(s in existingSamplers) {
            MCMCconf$addSampler(target = s$target,
                                type = "sampler_nonzero_indicator",
                                control = list(sampler = s,
                                               indicator = thisIndicatorNode),
                                silent = TRUE)
         }
      }
   }
   MCMCconf
}

