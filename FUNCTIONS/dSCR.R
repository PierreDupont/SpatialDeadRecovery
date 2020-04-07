#' @title Function to create a NIMBLE custom distribution with internalized detection
#'  probabilities calculation for faster SCR model runs.
#'
#' @description
#' \code{dSCR} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} of length n.detectors containing observation/non-observations 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors*2 with detectors x and y coordinates.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dSCR(p0 , sigma, sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], z[i,t] == 2)

#### 1.Density function ####
dSCR <- nimbleFunction(run = function( x = double(1)
                                     , pZero = double(0)
                                     , sigma = double(0)
                                     , sxy = double(1)
                                     , detector.xy = double(2)
                                     , indicator = double(0, default = 1.0)
                                     , log = integer(0, default = 0)){
   # Return type declaration
   returnType(double(0))
   outProb <- 0.0
   
   ## Check input dimensions
   n.detectors <- length(x)
   if(n.detectors <= 0){stop("invalid number of detectors")}
   numCoords <- length(sxy)
   if(numCoords <= 0){stop("invalid number of coordinates")}
   if(pZero < 0.0 | pZero > 1.0){outProb <- -Inf}
   if(sigma <= 0.0){outProb <- -Inf}
   if(dim(detector.xy)[1] != n.detectors){stop("detector coordinates has invalid dimension structure: invalid number of detectors")}
   if(dim(detector.xy)[2] != numCoords){stop("detector coordinates has invalid dimension structure: invalid number of coordinates")}
   
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
   
   ## Calculate the distances between AC and detectors
   d2 <- calculateDistance(sxy, detector.xy)
   
   ## Calculate the likelihood of the detection observations
   alpha <- -1.0 / (2.0 * sigma * sigma)
   detectLogLikeli <- rep(-Inf, n.detectors)
   for(j in 1:n.detectors){
      # Calculate the detection probability (negative distances repesent zero detection probability)
      p <- pZero * exp(alpha * d2[j])
      # Calculate the log-likelihood of each detection observation
      if(x[j] == 0){detectLogLikeli[j] <- log(1.0 - p)} else {detectLogLikeli[j] <- log(p)}
      # If the probability of detecting the current cell is zero then stop calculating and return the zero probability
      if(detectLogLikeli[j] <= -Inf){if(log == 0){return(0.0)}else{return(-Inf)}}
      }#j
   
   # Output
   outProb <- sum(detectLogLikeli[1:n.detectors])
   if(log == 0){outProb <- exp(outProb)}
   return(outProb)
   })

#### 2.Sampling function ####
rSCR <- nimbleFunction(run = function( n = integer(0)
                                     , pZero = double(0)
                                     , sigma = double(0)
                                     , sxy = double(1)
                                     , detector.xy = double(2)
                                     , indicator = double(0, default = 1.0)){
   # Return type declaration
   returnType(double(1))
   
   # Check input dimensions
   if(n!=0){print("rinhomPP only allows n = 1; using n = 1")}
   n.detectors <- dim(detector.xy)[1]
   if(n.detectors <= 0){stop("invalid number of detectors")}
   numCoords <- length(sxy)
   if(numCoords <= 0){stop("invalid number of coordinates")}
   if(pZero < 0.0 | pZero > 1.0){stop("invalid value for p0")}
   if(sigma <= 0.0){stop("invalid value for sigma")}
   if(dim(detector.xy)[2] != numCoords){stop("matrix of detector coordinates has invalid dimension structure")}
   
   ## Shortcut if individual is not available for detection
   if(indicator == 0){return(rep(0.0, n.detectors))}
   
   ## Calculate the distances between ACs and detectors
   D2 <- calculateDistance(sxy, detector.xy)
   
   ## Simulate the detections using the calculated detection distances ----
   alpha <- -1.0 / (2.0 * sigma * sigma)
   # Initialise a detections output vector
   detectOut <- rep(0, n.detectors)
   for(j in 1:n.detectors){
      p <- pZero * exp(alpha * D2[j])
      # Draw from a Bernoulli distribution with the calculated probability
      detectOut[j] <- rbinom(1, 1, p)
      }#j
   
   # Output
   return(detectOut)
})

#### 3.Registration ####
registerDistributions(list(
   dSCR = list(
      BUGSdist = "dSCR(pZero, sigma, sxy, detector.xy, indicator)",
      types = c( "value = double(1)", "pZero = double(0)", "sigma = double(0)"
                 ,"sxy = double(1)", "detector.xy = double(2)", "indicator = double(0)"),
      pqAvail = FALSE)))

