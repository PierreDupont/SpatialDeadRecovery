#' @title Function to create an Inhomogenous Poisson Process NIMBLE custom distribution.
#'
#' @description
#' \code{dinhomPP} Inhomogeneous Poisson Process to calculate the probability of a given AC location given an intensity surface.
#' 
#' @param x \code{Vector} of length 2. Coordinate values to calculate the density
#' @param intensityValues \code{Vector} of length n.cells denoting Values used in the intensity surface
#' @param intensityDims A \code{Vector} of length 2 denoting the dimensionality of the intensity surface
#' @param lowerCoords A \code{Vector} of length 2 denoting the lower coordinate values of the intensity surface
#' @param upperCoords A \code{Vector} of length 2 denoting the upper coordinate values of the intensity surface
#' @param log A \code{integer} required argument. It will always be log = TRUE when called from a model.
#'
#' @examples
#' sxy[i,1:2,t] ~ dinhomPP(rep(1,100), c(10,10), c(0,0), c(10,10)) 

#### 1.Density function ####
dinhomPP <- nimbleFunction(run = function( x = double(1)
                                         , intensityValues = double(1)
                                         , intensityDims = double(1)
                                         , lowerCoords = double(1)
                                         , upperCoords = double(1)
                                         , log = integer(0, default = 0)){
   
   # Return type declaration
   returnType(double(0))
   # Assess the dimensionality of the input coordinates
   dimCoords <- length(x)
   # Ensure that the dimensionality is valid
   if(dimCoords <= 0){stop("invalid dimension structure for the input coordinates")}
   # Recycle the input elements so that they are the same length as the dimensional coordinate values
   inIntensityDims <- numeric(length = dimCoords, value = intensityDims, recycle = TRUE)
   inLowerCoords <- numeric(length = dimCoords, value = lowerCoords, recycle = TRUE)
   inUpperCoords <- numeric(length = dimCoords, value = upperCoords, recycle = TRUE)
   inIntensityValues <- numeric(length = prod(inIntensityDims), value = intensityValues, recycle = TRUE)
   # Ensure that upper and lower coordinates are correctly oriented
   if(sum(inLowerCoords >= inUpperCoords) > 0) {stop("upper coordinates must be greater than the lower coordinates")}
      
   ## Calculate the dimension index offset
   # R's cumsum function is not implemented in nimbleFunction run-time code so we manually roll it out here using a for loop
   dimOffset <- integer(length = dimCoords)
   dimOffset[1] <- 1
   if(dimCoords > 1) {
      for(dimIter in 2:dimCoords) {
         # Calculate the dimensional index offset for each dimension
         dimOffset[dimIter] <- dimOffset[dimIter - 1] * inIntensityDims[dimIter - 1]
         }
      }
   ## Calculate the cell coordinates and indeces 
   # Calculate the zero-indexed cell coordinates of the input coordinates
   cellCoordinates <- trunc(floor((x - inLowerCoords) * inIntensityDims / (inUpperCoords - inLowerCoords)))
   outProb <- -Inf
   # Ensure that the calculated coordinates fall within the bounds of the intensity surface
   if(all(cellCoordinates >= 0 & cellCoordinates < inIntensityDims)) {
      # Calculate the cell index corresponding to the input coordinates
      cellIndex <- sum(cellCoordinates * dimOffset) + 1
      # Calculate the log probability density
      outProb <- log(inIntensityValues[cellIndex]) - log(sum(inIntensityValues)) +     # Relative probability of cell membership
      log(length(inIntensityValues)) - sum(log(inUpperCoords - inLowerCoords))      # Normalisation by cell area
      }
   
   ## Return the retrieved density ----
   if(log == 0){outProb <- exp(outProb)}
   return(outProb)
   })

#### 2.Sampling function ####
rinhomPP <- nimbleFunction(run = function( n = integer(0)
                                         , intensityValues = double(1)
                                         , intensityDims = double(1)
                                         , lowerCoords = double(1)
                                         , upperCoords = double(1)){
   

   # Return type declaration
   returnType(double(1))
   
   ## Sanity test the inputs ----
   # Assess the dimensionality of the input coordinates
   dimCoords <- length(intensityDims)
   # Ensure that the dimensionality is valid
   if(dimCoords <= 0) {
      stop("invalid dimension structure for the input coordinates")
      }
   # Recycle the input elements so that they are the same length as the dimensional coordinate values
   inIntensityDims <- numeric(length = dimCoords, value = intensityDims, recycle = TRUE)
   inLowerCoords <- numeric(length = dimCoords, value = lowerCoords, recycle = TRUE)
   inUpperCoords <- numeric(length = dimCoords, value = upperCoords, recycle = TRUE)
   inIntensityValues <- numeric(length = prod(inIntensityDims), value = intensityValues, recycle = TRUE)
   # Ensure that upper and lower coordinates are correctly oriented
   if(sum(inLowerCoords >= inUpperCoords) > 0){stop("upper coordinates must be greater than the lower coordinates")}
   # Ensure that only one sample is requested
   if(n <= 0) {
      stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rinhomPP only allows n = 1; using n = 1")
      }
   
   ## Calculate the dimension index offset 
   # R's cumsum function is not implemented in nimbleFunction run-time code so we manually roll it out here using a for loop
   dimOffset <- numeric(length = dimCoords)
   dimOffset[1] <- 1
   if(dimCoords > 1) {
      for(dimIter in 2:dimCoords) {
         # Calculate the dimensional index offset for each dimension
         dimOffset[dimIter] <- dimOffset[dimIter - 1] * inIntensityDims[dimIter - 1]
         }
      }
   
   ## Sample a random set of cell coordinates from the intensity values 
   # Sample the cell index directly
   cellIndex <- rcat(1, prob = inIntensityValues)
   # Convert the cell index into a set of multidimensional cell coordinates (zero indexed)
   cellCoordinates_uncycled <- trunc(floor(rep(cellIndex - 1, dimCoords) / dimOffset))
   cellCoordinates <- numeric(length = dimCoords)
   for(coordIter in 1:dimCoords) {
      # According to the Nimble Manual (version 0.6-8) the modulo operator has not been defined to accept
      # vector input in nimbleFunction run-time code so, for safety, it has been rolled out into a 'for'
      # loop here.  It may be possible in furture versions of Nimble to simply replace the loop with
      # 'cellCoordinates <- cellCoordinates_uncycled %% intensityDims'.
      cellCoordinates[coordIter] <- cellCoordinates_uncycled[coordIter] %% inIntensityDims[coordIter]
      }
   
   ## Draw a random coordinate within the cell coordinates ----
   # Calculate the cell size in each dimension
   cellsize <- (inUpperCoords - inLowerCoords) / inIntensityDims
   # Retrieve the output coordinates
   outCoordinates <- inLowerCoords + (cellCoordinates + runif(dimCoords, min = 0.0, max = 1.0)) * cellsize
   return(outCoordinates)
   })

#### 3.Registration ####
registerDistributions(list(
   dinhomPP = list( BUGSdist = "dinhomPP(intensityValues, intensityDims, lowerCoords, upperCoords)"
                  , types = c("value = double(1)", "intensityValues = double(1)", "intensityDims = double(1)", "lowerCoords = double(1)", "upperCoords = double(1)")
                  , pqAvail = FALSE)))

