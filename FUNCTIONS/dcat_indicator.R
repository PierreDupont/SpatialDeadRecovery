#' @title Function to create a NIMBLE categorical distribution with indicator for faster SCR model runs.
#'
#' @description
#' \code{dcat_indicator} returns the likelihood of being observed in a given category/state/trap/detector
#' 
#' @param x \code{numeric} denoting the observed category (0 stands for unobserved = not sampled).
#' @param p A \code{vector} of length n.categories denoting the probabilities to belong to each category.
#' @param indicator A \code{numeric}; if == 0 no observation possible ; category == 0 ; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' dead.recovery[i,t] ~ dcat_indicator(p[1:n.detectors], z[i,t] == 2)

#### 1.Density function ####
dcat_indicator <- nimbleFunction(run = function( x = double(0)
                                               , p = double(1)
                                               , indicator = double(0, default = 1.0)
                                               , log = integer(0, default = 0)){
   # Return type declaration
   returnType(double(0))
   ## Shortcut if individual is not available for detection
   if(indicator == 0){
      if(x == 0){
         if(log == 0) return(1.0) else return(0.0)
      }else{
         if(log == 0) return(0.0) else return(-Inf)
      }
   }

   # Calculate the log-likelihood of category x over n.categories 
   logProb <- dcat(x, prob = p, log = TRUE)
   if(log)return(logProb)
   return(exp(logProb))
   })

#### 2.Sampling function ####
rcat_indicator <- nimbleFunction(run = function( n = integer(0)
                                               , p = double(1)
                                               , indicator = double(0, default = 1.0)){
   # Return type declaration
   returnType(double(0))
   # Check input dimensions
   if(n!=1){print("rinhomPP only allows n = 1; using n = 1")}
   n.categories <- length(p)
   ## Shortcut if individual is not available for detection
   if(indicator == 0){return(0.0)}
   # Draw from a Categorical distribution 
   return(rcat(1, p[1:n.categories]))
   })

#### 3.Registration ####
registerDistributions(list(
   dcat_indicator = list(
      BUGSdist = "dcat_indicator(p, indicator)",
      types = c( "value = double(0)", "p = double(1)", "indicator = double(0)"),
      pqAvail = FALSE)))

