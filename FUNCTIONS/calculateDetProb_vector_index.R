#' @title NIMBLE Function to calculate the sqaure distance between an AC location and several detectors.
#'
#' @description
#' \code{calculateDetProb_vector} is a NIMBLE Function to calculate the detector-specific detection probability
#' following a half-normal detection function of the distance between the AC location and detectors.
#' 
#' @param sxy \code{vector} of length 2 containing x and y coordinates of the individual activity center. 
#' @param detector.xy \code{matrix} of dimensions n.detectors*2 denoting the x and y detectors coordinates.
#' @param pZero \code{vector} of length n.detectors denoting the detector-specific baseline detection probability.
#' @param sigma \code{numeric} denoting the scale parameter of the detection function (related to home-range size).
#' @param maxDist \code{numeric} denoting the maximum distance to be considered for potential detections (optional).
#'
#' @examples
#' p[i,1:n.detectors,t] <- calculateDetProb_vector(sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], p0[i,1:n.detectors,t], sigma[i,t], maxDist)

calculateDetProb_vector_index <- nimbleFunction(run = function( sxy = double(1)
                                                          , detector.xy = double(2)
                                                          , pZero = double(1)
                                                          , pZeroIndex = double(1)
                                                          , sigma = double(0)
                                                          , maxDist= double(0, default = 0.0)
                                                          , indicator = double(0, default = 1.0)){
  # Return type declaration
  returnType(double(1))
  
  # Check input dimensions
  n.detectors <- dim(detector.xy)[1]
  alpha <- -1 / (2 * sigma * sigma)
  
  p <- rep(0.0, n.detectors)
  if(indicator == 0){return(p)}
  
  # Calculate distance vector (Output)
  d2 <- pow(detector.xy[1:n.detectors,1] - sxy[1], 2) + pow(detector.xy[1:n.detectors,2] - sxy[2], 2)
  
  # Calculate detector-specific probabilities
  if(maxDist > 0){
    for(j in 1:n.detectors){
      if(d2[j] <= (maxDist*maxDist)){p[j] <- pZero[pZeroIndex[j]] * exp(alpha * d2[j])}
    }#j
  }else{p <- pZero * exp(alpha * d2)}
  
  return(p)
})
