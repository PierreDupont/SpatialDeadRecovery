#' @title Function to create a NIMBLE custom distribution for faster SCR model runs. i
#'
#' @description
#' \code{dSCR_LESS} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} of length n.detectors containing observation/non-observations 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors*2 with detectors x and y coordinates.
#' @param maxDist A \code{Numeric} with the maximum distance to the AC where detections are considered possible.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dSCR_LESS(p0 , sigma, sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], 10, z[i,t] == 2)

#### 1.Density function ####
dbin_LESSWolf <- nimbleFunction(run = function( x = double(1)
                                            , sigma = double(0)
                                            , trials = double(1)
                                            , d2 = double(1)
                                            , maxDist = double(0, default = 0.0)
                                            , indicator = double(0, default = 1.0)
                                            , z = double(0)
                                            , p0State = double(2)
                                            , idResponse = double(0)
                                            , detCountries = double(1)
                                            , detRoads = double(1)
                                            , detTracks = double(1)
                                            , betaResponse = double(0)  
                                            , betaTracks = double(0)
                                            , betaRoads = double(0)
                                            , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  ## Check input dimensions
  n.detectors <- length(x)
  pZero <- numeric(n.detectors, init = FALSE)
  
  ## Shortcut if individual is not available for detection
  if(indicator == 0){
    if(sum(x[1:n.detectors]) == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  
  
  for(j in 1:n.detectors){
      pZero[j] <- ilogit(logit(p0State[detCountries[j], z]) +
                                     betaResponse * idResponse +
                                     betaRoads * detRoads[j] +
                                     # betaSnow * detSnow[j,t] +
                                     betaTracks * detTracks[j])
  }#j
  
  
  
  
  ## Calculate the likelihood of the detection observations
  alpha <- -1.0 / (2.0 * sigma * sigma)
  detectLogLikeli <- numeric(n.detectors, init = FALSE)
  p <- numeric(n.detectors, init = FALSE)
  
  maxDist_squared <- maxDist*maxDist
  for(j in 1:n.detectors){
    # Calculate the detection probability (negative distances repesent zero detection probability)
    if(d2[j] <= (maxDist_squared)){p[j] <- pZero[j] * exp(alpha * d2[j])} else {p[j] <- 0.0000000}
    
    # # Calculate the log-likelihood of each detection observation
    # if(x[j] == 0){detectLogLikeli[j] <- log(1.0 - p)} else {detectLogLikeli[j] <- log(p)}
    # # If the probability of detecting the current cell is zero then stop calculating and return the zero probability
    # if(detectLogLikeli[j] <= -Inf){if(log == 0){return(0.0)}else{return(-Inf)}}
  }#j
  
  
  logProb <- sum(dbinom(x, prob = p, size = trials, log = TRUE))
  if(log)return(logProb)
  return(exp(logProb))
  
  
  # Output
  # outProb <- sum(detectLogLikeli[1:n.detectors])
  # if(log == 0){outProb <- exp(outProb)}
  # return(outProb)
})

#### 2.Sampling function ####
rbin_LESSWolf <- nimbleFunction(run = function( n = integer(0)
                                            , sigma = double(0)
                                            , trials = double(1)
                                            , d2 = double(1)
                                            , maxDist = double(0, default = 0.0)
                                            , indicator = double(0, default = 1.0)
                                            , z = double(0)
                                            , p0State = double(2)
                                            , idResponse = double(0)
                                            , detCountries = double(1)
                                            , detRoads = double(1)
                                            , detTracks = double(1)
                                            , betaResponse = double(0)  
                                            , betaTracks = double(0)
                                            , betaRoads = double(0)){
  # Return type declaration
  returnType(double(1))
  
  # Check input dimensions
  n.detectors <- length(d2)
  
  ## Shortcut if individual is not available for detection
  if(indicator == 0){return(rep(0.0, n.detectors))}
  
  pZero <- numeric(n.detectors, init = FALSE)
  
  for(j in 1:n.detectors){
    pZero[j] <- ilogit(logit(p0State[detCountries[j], z]) +
                         betaResponse * idResponse +
                         betaRoads * detRoads[j] +
                         # betaSnow * detSnow[j,t] +
                         betaTracks * detTracks[j])
  }#j
  ## Simulate the detections using the calculated detection distances ----
  alpha <- -1.0 / (2.0 * sigma * sigma)
  # Initialise a detections output vector
  detectOut <- rep(0, n.detectors)
  maxDist_squared <- maxDist*maxDist
  for(j in 1:n.detectors){
    # Calculate the detection probability (negative distances repesent zero detection probability)
    if(d2[j] <= maxDist_squared){
      p <- pZero[j] * exp(alpha * d2[j])
      # Draw from a Bernoulli distribution with the calculated probability
      detectOut[j] <- rbinom(1, trials[j], p)
    }#if
  }#j
  
  # Output
  return(detectOut)
})

#### 3.Registration ####
registerDistributions(list(
  dbin_LESSWolf = list(
    BUGSdist = "dbin_LESSWolf(sigma, trials, d2, maxDist, indicator, z, p0State, idResponse,
    detCountries, detRoads, detTracks, betaResponse, betaTracks, betaRoads )",
    types = c( "value = double(1)", "sigma = double(0)", "trials = double(1)", "d2 = double(1)", 
               "maxDist = double(0)" ,"indicator = double(0)", "z = double(0)", "p0State = double(2)",
                 "idResponse = double(0)", "detCountries = double(1)", "detRoads = double(1)","detTracks = double(1)"
                 , "betaResponse = double(0)" , "betaTracks = double(0)", "betaRoads = double(0)"),
    pqAvail = FALSE)))

