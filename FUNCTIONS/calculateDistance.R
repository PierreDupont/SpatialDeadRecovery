#' @title NIMBLE Function to calculate the square distance between an AC location and a vector of detectors locations.
#'
#' @description
#' \code{calculateDistance} is a NIMBLE Function to calculate the square distance between an AC location and several detectors.
#' 
#' @param sxy \code{Vector} of length 2 containing location X and Y of the activity center of one individual . 
#' @param detector.xy \code{Matrix} with nrow = n.detectors and ncol the X and Y coordinates.
#' @param indicator \code{numeric} denotes whether the individual is considered alive (1) or dead (0) and whether the  square distance should be calculated or not. .
#'
#' @examples
#' d2[i,1:n.detectors,t] <- calculateDistance(sxy[i,1:2,t], detector.xy[1:n.detectors,1:2], z[i,t]==2)

calculateDistance <- nimbleFunction(run = function( sxy = double(1)
                                                  , detector.xy = double(2)
                                                  , indicator = double(0, default = 1.0)){
   # Return type declaration
   returnType(double(1))
      
   # Check input dimensions
   n.detectors <- dim(detector.xy)[1]
   if(indicator == 0){
      return(rep(0,n.detectors))
   }
   # Calculate distance vector (Output)
   d2 <- pow(detector.xy[1:n.detectors,1] - sxy[1], 2) + pow(detector.xy[1:n.detectors,2] - sxy[2], 2)
   # Return square distance
   return(d2)
   })
