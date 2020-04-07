#' @title Function to create a NIMBLE Vectorized binomial custom distribution.
#'
#' @description
#' \code{dbin_vector} Vectorized binomial distribution to calculate probability of all detections for one individual at once.
#' Note that if the indicator is 0, we can just return logProb = 0.  This corresponds to
#' a virtual individual who is not in the model, so probability of detection is 0 and the fake data will be all 0s.
#' The model declaration  y[i, 1:n.detectors] ~ dbin_vector(p[i, 1:n.detectors], z[i], trials)
#' will become a call effectively like logProb_y[i] <- dbin_vector(y[i, 1:n.detectors], p[i, 1:n.detectors], z[i], trials)
#' 
#' @param x \code{Vector} of length n.detectors containing observation/non-observations. 
#' @param prob \code{Vector} of length n.detectors denoting the detector-specific detection probability.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param trials A \code{Numeric} denoting the number of trials for each observation (the same for all observations so far)
#' @param log A \code{integer} required argument.  It will always be log = TRUE when called from a model.
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dbin_vector(p0 , sigma, sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], z[i,t] == 2)

#### 1.Density function ####
dbin_vector <- nimbleFunction( run = function( x = double(1)
                                             , prob = double(1)
                                             , indicator = double(0, default = 1.0)
                                             , trials = double(1)
                                             , log = integer(0, default = 0)){
   returnType(double(0))
   if(indicator == 0){
      if(sum(x) == 0){
         if(log == 0) return(1.0) else return(0.0)
         }else{ if(log == 0) return(0.0) else return(-Inf) }
      }#if
   logProb <- sum(dbinom(x, prob = prob, size = trials, log = TRUE))
   if(log) return(logProb)
   return(exp(logProb))
   })

#### 2.Sampling function ####
rbin_vector <- nimbleFunction(run = function( n = double(0, default = 1.0)
                                            , prob = double(1)
                                            , indicator = double(0, default = 1.0)
                                            , trials = double(1)){
   returnType(double(1))
   n.detectors <- length(prob)
   if(indicator == 0){return(rep(0.0, n.detectors))}
   detectOut <- rbinom(n = n.detectors, size = trials, prob = prob)
   return(detectOut)
   })

#### 3.Registration ####
registerDistributions(list(
   dbin_vector = list(
      BUGSdist = "dbin_vector(prob, indicator, trials)",
      types = c("value = double(1)", "prob = double(1)", "indicator = double(0)", "trials = double(1)"),
      pqAvail = FALSE)))