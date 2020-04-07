#' @title NIMBLE Function to calculate the sqaure distance between an AC location and several detectors.
#'
#' @description
#' \code{calculateDistance} is a NIMBLE Function to calculate the sqaure distance between an AC location and several detectors.
#' 
#' @param sxy \code{Vector} of length 2 denoting an individual's activity center coordinates. 
#' @param detector.xy \code{Matrix} of dimension n.detectors*2 denoting the coordinates of detectors.
#' @param pZero \code{numeric} denoting the baseline detection probability.
#' @param sigma \code{numeric} denoting the scale parameter of the detection function.
#' @param maxDist \code{numeric} denoting the maximum distance where detections are allowed (saves time).
#' @param indicator \code{numeric} denoting whether the individual is available for detection or not (if not p = 0 ; saves time).
#'
#' @examples
#' p[i,1:n.detectors,t] <- calculateDetProb(sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], p0[i,t], sigma, 20, z[i,t] == 2)

calculateDetProbDead <- nimbleFunction(run = function( sxy = double(1)
                                                 , detector.xy = double(2)
                                                 , sigma = double(0)
                                                 , maxDist = double(0, default = 0.0)
                                                 , indicator = double(0, default = 1.0)){
   # Return type declaration
   returnType(double(1))
   
   # Check input dimensions
   n.detectors <- dim(detector.xy)[1]
   
   # Initialize p vector
   p <- rep(0.0, n.detectors)
   if(indicator == 0){return(p)}
   
   # Calculate distance vector (Output)
   alpha <- -1/(2 * sigma * sigma)
   d2 <- pow(detector.xy[1:n.detectors,1] - sxy[1], 2) + pow(detector.xy[1:n.detectors,2] - sxy[2], 2)
   
   # Calculate detector-specific probabilities
   if(maxDist > 0){
      for(j in 1:n.detectors){
         if(d2[j] <= (maxDist*maxDist)){p[j] <- exp(alpha * d2[j])}
         }#j
      }else{
      p <- exp(alpha * d2)
      }
   
   return(p)
   })
