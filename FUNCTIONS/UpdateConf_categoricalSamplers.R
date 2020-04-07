#' @title NIMBLE script to create a new categorical NIMBLE sampler.
#'
#' @description
#' \code{sampler_categorical_faster} creates a faster categorical sampler.
#' It avoids calculations of categories where the prior probability is 0.
#'
#' \code{updateConf_DASamplers} modifies the current NIMBLE configuration to include faster categorical sampler. 
#'
#' @examples
#' MCMCconf <- updateCategoricalSamplers(MCMCconf = MCMCconf)

#### 1.Faster Categorical Sampler ####
sampler_categorical_faster <- nimbleFunction(
   name = 'sampler_categorical_faster',
   contains = sampler_BASE,
   setup = function(model, mvSaved, target, control) {
      ## node list generation
      targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
      calcNodesNoSelf  <- model$getDependencies(target, self = FALSE)
      ## numeric value generation
      k <- length(model$getParam(target, 'prob'))
      probs <- numeric(k)
      ## checks
      if(length(targetAsScalar) > 1)  stop('cannot use categorical sampler on more than one target node')
      if(model$getDistribution(target) != 'dcat') stop('can only use categorical sampler on node with dcat distribution')
   },
   run = function() {
      currentValue <- model[[target]]
      probs[currentValue] <<- exp(getLogProb(model, target) + getLogProb(model, calcNodesNoSelf))
      for(i in 1:k) {
         if(i != currentValue) {
            model[[target]] <<- i
            priorLogProb <- calculate(model, target)
            if(priorLogProb == -Inf | is.nan(priorLogProb))
               probs[i] <<- 0
            else {
               probs[i] <<- exp(priorLogProb + calculate(model, calcNodesNoSelf))
               if(is.nan(probs[i])) probs[i] <<- 0
            }
         }
      }
      newValue <- rcat(1, probs)
      if(newValue != currentValue) {
         model[[target]] <<- newValue
         calculate(model, target)
         calculate(model, calcNodesNoSelf)
         nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
         nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelf, logProb = TRUE)
      } else {
         nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
         nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelf, logProb = TRUE)
      }
   },
   methods = list(
      reset = function() { }
   )
)

#### 2.Update function ####
updateCategoricalSamplers <- function(MCMCconf) {
   samplers <- MCMCconf$getSamplers()
   samplerNames <- unlist(lapply(samplers, `[[`, 'name'))
   boolCategorical <- samplerNames == 'categorical'
   targets <- unlist(lapply(samplers, `[[`, 'target'))
   categoricalTargets <- targets[boolCategorical]
   MCMCconf$removeSamplers(categoricalTargets)
   for(target in categoricalTargets) {
      MCMCconf$addSampler(target = target,
                          type = "sampler_categorical_faster")
   }
   MCMCconf
}