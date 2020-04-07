#' @title NIMBLE function to check whether the sxy coordinates falls within the spatial domain and within suitable habitat
#'
#' @description
#' \code{OKHab} returns the 1 if the sxy coordinates falls within the spatial domain and within suitable habitat, 0 otherwise
#' 
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param lowerCoords A \code{Vector} of length 2 with lower x and y coordinates of the spatial domain.
#' @param upperCoords A \code{Vector} of length 2 with upper x and y coordinates of the spatial domain.
#' @param habitatmx A \code{Matrix} with representing habitat suitability
#'
#' @examples
#' pOK[i,t] <- OKHab(sxy[i,1:2,t], lowerCoords[1:2], upperCoords[1:2], habitat.mx[1:y.max,1:x.max])

#### 1.Density function ####
OKHab <- nimbleFunction(run = function(    sxy = double(1)
                                         , lowerCoords = double(1)
                                         , upperCoords = double(1)
                                         , habitatmx = double(2)){
  # Return type declaration
  returnType(double(0))
  
  ## Check input dimensions
  
  
  ## If sxy outside of the habitat min/max
  if(((sxy[1] >  lowerCoords[1])  *  (sxy[1] <  upperCoords[1]) * 
      (sxy[2] >  lowerCoords[2])  *  (sxy[2] <  upperCoords[2]))==0 ){
    return(0.0) 
  }else{
    pOK <- habitatmx[trunc(sxy[2])+1, trunc(sxy[1])+1]
    
     return(pOK) 
  }
  
})