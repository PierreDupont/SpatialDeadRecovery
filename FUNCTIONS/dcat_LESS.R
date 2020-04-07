#' @title Function to create a NIMBLE custom distribution for faster SCR model runs. 
#'
#' @description
#' \code{dcat_LESS} returns the likelihood of a given categorical spatial detection y[i] 
#'
#' @param x \code{Vector} of length n.detectors containing observations/non-observations (dbern)
#' @param pZero \code{numeric} denoting the baseline detection probability of the half-normal.
#' @param sigma \code{Numeric} scalar denoting the scale parameter of the half-normal.
#' @param d2 \code{Vector} denoting the distance to each detector.
#' @param maxDist \code{Numeric} with the maximum distance to the AC where detections are considered possible.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{Integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,t] ~ dcat_LESS_vector(p0[1:n.detectors], sigma, d2[i,1:n.detectors,t], 10, z[i,t]==3)

#### 1.Density functions ####
dcat_LESS <- nimbleFunction(run = function( x = double(0),
                                               pZero = double(0),
                                               sigma = double(0),
                                               d2 = double(1),
                                               maxDist = double(0, default = 0.0),
                                               indicator = double(0, default = 1.0),
                                               log = integer(0, default = 0)){
   ## Return type declaration
   returnType(double(0))
   
   ## Check input dimensions
   n.detectors <- length(d2)
   
   ## Shortcut if individual is not available for detection
   if(indicator == 0){
      if(x == 0){
         if(log == 0) return(1.0)
         else return(0.0)
      } else {
         if(log == 0) return(0.0)
         else return(-Inf)
      }
   }
   
   ## Process input values
   alpha <- -1.0 / (2.0 * sigma * sigma)
   maxDist_squared <- maxDist*maxDist
   
   ## Calculate detector-specific detection probabilities (only within maxDist of the AC ; else p = 0)
   p <- numeric(length = n.detectors, value = 0)
   for(j in 1:n.detectors){
      if(d2[j] <= maxDist_squared){p[j] <- pZero * exp(alpha * d2[j])}
   }#j
   
   # Check that there is at least one possible detector with p > 0.0
   if(sum(p[1:n.detectors]) <= 0.0){      
      if(log == 0) return(0.0)
      else return(-Inf)
   }
   
   # Calculate the log-likelihood of category x over n.categories 
   logProb <- dcat(x, prob = p, log = TRUE)
   if(log)return(logProb)
   return(exp(logProb))
})

#### 2.Sampling functions ####
rcat_LESS <- nimbleFunction(run = function( n = integer(0),
                                            pZero = double(0),
                                            sigma = double(0),
                                            d2 = double(1),
                                            maxDist = double(0, default = 0.0),
                                            indicator = double(0, default = 1.0)){
   ## Return type declaration
   returnType(double(0))
   
   ## Check input dimensions
   n.detectors <- length(d2)
   
   ## Shortcut if individual is not available for detection
   if(indicator == 0){return(0.0)}
   
   ## Process input values
   alpha <- -1.0 / (2.0 * sigma * sigma)
   maxDist_squared <- maxDist * maxDist
   
   ## Calculate detector-specific detection probabilities (only within maxDist of the AC ; else p = 0)
   p <- numeric(length = n.detectors, value = 0)
   for(j in 1:n.detectors){
      if(d2[j] <= maxDist_squared){p[j] <- pZero * exp(alpha * d2[j])}
   }#j
   
   ## Draw detections from a Categorical distribution with the calculated probabilities
   y <- rcat(n = 1, prob = p)
   
   ## Output
   return(y)
})

#### 3.Registration ####
registerDistributions(list(
   dcat_LESS = list(
      BUGSdist = "dcat_LESS(pZero, sigma, d2, maxDist, indicator)",
      types = c( "value = double(0)", "pZero = double(0)", "sigma = double(0)",
                 "d2 = double(1)", "maxDist = double(0)", "indicator = double(0)"),
      pqAvail = FALSE)))