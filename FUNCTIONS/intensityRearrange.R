#' @title NIMBLE function to rearrange the intensity surface for use in NIMBLE
#'
#' @description
#' \code{intensityRearrange} is a NIMBLE function to rearrange the intensity surface for use in NIMBLE.
#' 
#' @param intensityValues \code{Vector} of length n.cells containing the inensity surface values. 
#' @param intensityDims \code{Vector} of length 2 denoting the dimenions of the space considered.
#' @param dimOrder \code{Vector} of length 2 denoting ????.
#' @param dimBackwards \code{Vector} of length 2 denoting ????.
#'
#' @examples
#' mu[1:n.cells] <- intensityRearrange(OldMu[1:n.cells], c(0,10), c(0,0), c(0,0))

intensityRearrange <- nimbleFunction(run = function( intensityValues = double(1)
                                                   , intensityDims = double(1)
                                                   , dimOrder = double(1, default = c(0))
                                                   , dimBackwards = integer(1, default = c(0))){
   ## Specify the return type dimensionality 
   returnType(double(1))
   ## Sanity test the inputs 
   # Assess the dimensionality of the input coordinates
   dimCoords <- length(intensityDims)
   # Ensure that the dimensionality is valid
   if(dimCoords <= 0){stop("invalid dimension structure for the input coordinates")}
   # Recycle the input element so that they are the same length as the dimensional coordinate values
   inIntensityDims <- numeric(length = dimCoords, value = intensityDims, recycle = TRUE)
   numCells <- prod(inIntensityDims)
   inDimOrder <- numeric(length = dimCoords, value = dimOrder, recycle = TRUE)
   inDimBackwards <- integer(length = dimCoords, value = dimBackwards, recycle = TRUE)
   inIntensityValues <- numeric(length = numCells, value = intensityValues, recycle = TRUE)
   ## Calculate the dimension index offset (new dimension structure) 
   # R's cumsum function is not implemented in nimbleFunction run-time code so we manually roll it out here using a for loop
   dimOffsetNew <- numeric(length = dimCoords)
   dimOffsetNew[1] <- 1
   if(dimCoords > 1) {
      for(dimIter in 2:dimCoords){
         dimOffsetNew[dimIter] <- dimOffsetNew[dimIter - 1] * inIntensityDims[dimIter - 1]
         }#dimIter
      }#if
   ## Calculate the dimension index offset (original dimension structure) 
   # According to NIMBLE's manual (version 0.6-8) the 'sort', 'rank', 'ranked' and 'order' functions have not yet
   # been implemented in run code.  I've implemented a clunky sort here.  This won't be very efficient for
   # high-dimensional problems.
   dimRank <- integer(length = dimCoords, value = rep(0, dimCoords))
   curVal <- min(inDimOrder)
   curRank <- 1
   # Create a vector of order integers in dimRank by ranking the values in inDimOrder.  Ties are
   # seperated by the order that they occur in the dimensional array.
   while(curRank <= dimCoords) {
      minVal <- curVal
      for(dimIter in 1:dimCoords) {
         if(inDimOrder[dimIter] == curVal) {
            # Set the dimension to the current rank
            dimRank[dimIter] <- curRank
            # Increment the rank
            curRank <- curRank + 1
         } else if(inDimOrder[dimIter] > curVal) {
            # Otherwise find the next current 
            if(minVal == curVal) {
               minVal <- inDimOrder[dimIter]
            } else {
               minVal <- min(minVal, inDimOrder[dimIter])
            }
         }
      }
      curVal <- minVal
   }
   # Iterate over the ranks and calculate the dimensional offset of each dimension corresponding to the rank
   dimOffsetOld <- numeric(length = dimCoords)
   curDimOffset <- 1
   for(rankIter in 1:dimCoords) {
      for(dimIter in 1:dimCoords) {
         if(dimRank[dimIter] == rankIter) {
            dimOffsetOld[dimIter] <- curDimOffset
            curDimOffset <- curDimOffset * inIntensityDims[dimIter]
            }
         }
      }
   ## Convert each old array index to a new array index
   # Create a vector to store the cell coordinates
   curCellCoords <- numeric(length = dimCoords)
   # Create a vector to store the rearranged intensity values
   newIntensity <- numeric(length = numCells)
   for(oldIndex in 1:numCells) {
      # Convert the cell index into a set of multidimensional cell coordinates (zero indexed)
      curCellCoords_uncycled <- trunc(floor(rep(oldIndex - 1, dimCoords) / dimOffsetOld))
      for(coordIter in 1:dimCoords) {
         curCellCoords[coordIter] <- curCellCoords_uncycled[coordIter] %% inIntensityDims[coordIter]
         if(inDimBackwards[coordIter] != 0) {
            # Reverse the dimensional coordinate of the cell if the coordinate system is reversed in this dimension
            curCellCoords[coordIter] <- inIntensityDims[coordIter] - 1 - curCellCoords[coordIter]
            }
         }
      # Create the new cell index from the coordinates
      newIndex <- sum(curCellCoords * dimOffsetNew) + 1
      # Copy the element from the old intensity values vector to the new intensity values vector
      newIntensity[newIndex] <- inIntensityValues[oldIndex]
      }
   return(newIntensity)
   })