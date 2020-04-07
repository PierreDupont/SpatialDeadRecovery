#' @title NIMBLE script to create a reversible-jump NIMBLE sampler.
#'
#' @description
#' \code{sampler_DAindicator} creates a sampler that proposes to insert or remove entities labeled 
#' by an indicator node, along with associated variables that are sampled from their prior when an insertion proposal is made.
#' 
#' \code{updateConf_DASamplers} modifies the current NIMBLE configuration to include reversible jump insertion/removal 
#' of entities indicated by indicator variable and related variables.
#' The related variables will be sampled according to their priors when proposed for insertion.
#' 
#' @examples
#' MCMCconf <- updateConf_DASamplers(MCMCconf = MCMCconf, indicatorVar = "z", relatedVarList  = list("sxy", "p0", "sigma"))

#### 1.Reversible-Jump Sampler ####
sampler_DAindicator <- nimbleFunction( contains = sampler_BASE
                                       , setup = function(model, mvSaved, target, control){
                                          indicatorNode <- target
                                          relatedNodes <- control$related
                                          insertProb <- if(is.null(control$insertProb)) 0.2 else control$insertProb
                                          removeProb <- if(is.null(control$removeProb)) 0.2 else control$removeProb
                                          logInsertProb <- log(insertProb)
                                          logRemoveProb <- log(removeProb)
                                          calcNodes <- model$getDependencies(c(indicatorNode, relatedNodes))
                                          calcNodesReduced <- model$getDependencies(indicatorNode)
                                       }
                                       , run = function(){
                                          currentIndicator <- model[[indicatorNode]]
                                          currentLogProb <- getLogProb(model, calcNodes)
                                          if(currentIndicator == 1) {
                                             if(runif(1,0,1) > removeProb) return()
                                             ## propose removal
                                             currentLogProb <- model$getLogProb(calcNodes)
                                             logProbReverseProposal <- model$getLogProb(relatedNodes)
                                             model[[indicatorNode]] <<- 0
                                             calculate(model, calcNodes)
                                             log_accept_prob <-
                                                model$getLogProb(calcNodesReduced) -
                                                currentLogProb +
                                                logProbReverseProposal +
                                                logInsertProb - logRemoveProb
                                          } else {
                                             if(runif(1,0,1) > insertProb) return()
                                             ## propose insertion
                                             ## 1. Treat the relatedNodes as if mixing according to their prior
                                             currentLogProb <- model$getLogProb(calcNodesReduced)
                                             model$simulate(relatedNodes)
                                             model[[indicatorNode]] <<- 1
                                             ## jumping density
                                             proposalLogProb <- model$calculate(calcNodes)
                                             logProbForwardProposal <- model$getLogProb(relatedNodes)
                                             log_accept_prob <-
                                                proposalLogProb -
                                                currentLogProb -
                                                logProbForwardProposal +
                                                logRemoveProb - logInsertProb
                                          }
                                          jump <- decide(log_accept_prob)
                                          if(jump)
                                             copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
                                          else
                                             copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
                                       }
                                       , methods = list(reset = function(){}))

#### 2.Update function ####
updateConf_DASamplers <- function( MCMCconf
                                 , indicatorVar
                                 , relatedVarList
                                 , byRow = NULL
                                 , control = list()){
   model <- MCMCconf$model
   indicatorNodes <- model$expandNodeNames(indicatorVar)
   groups <- vector(length(indicatorNodes), mode = "list")
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
         newGroups <- matrix(relatedNodes, nrow = relatedDims[1], byrow = !byRow)
         newGroups <- apply(newGroups, 1, list)
      }
      if(relatedNDim == 1) {
         newGroups <- lapply(relatedNodes, list)
      }
      groups <- mapply(function(a, b) list(c(a[[1]], b[[1]])),
                       groups,
                       newGroups,
                       SIMPLIFY = FALSE)
   }
   boolData <- model$isData(indicatorNodes)
   indicatorNodes <- indicatorNodes[!boolData]
   groups <- groups[!boolData]
   for(i in seq_along(groups)) {
      thisIndicatorNode <- indicatorNodes[i]
      thisRelatedNodes <- groups[[i]][[1]]
      MCMCconf$removeSamplers(thisIndicatorNode)        
      MCMCconf$addSampler(target = thisIndicatorNode,
                          type = "sampler_DAindicator",
                          control = c(list(relatedNodes = thisRelatedNodes),
                                      control),
                          silent = TRUE)
   }
   MCMCconf
}