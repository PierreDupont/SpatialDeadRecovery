getDensity <- nimbleFunction(run = function( curCoords = double(2),
                                             lowerCoords = double(2),
                                             upperCoords = double(2),
                                             indicator = double(1)){
   
   ## 1. Specify the return type dimensionality ----
   # Return type declaration
   returnType(double(1))
   
   ## 2. Assess dimensions ----
   # Assess the dimensionality of the input coordinates
   dimCoords <- dim(lowerCoords)[2]
   # Assess the number of observation windows
   numWindows <- dim(lowerCoords)[1]
   # Assess the number of points
   numPoints <- dim(curCoords)[1]
   
   ## 3. Assess intersection between points and observation windows ----
   density <- numeric(length = numWindows, value = 0.0, recycle = TRUE)
   # The observation window is an area/volume so assess that it falls within that volume
   for(pointIter in 1:numPoints){
      if(indicator[pointIter] == 1){
         isInWindow <- numeric(length = numWindows, value = 1.0, recycle = TRUE)
         for(dimIter in 1:dimCoords){
            isInWindow <- isInWindow * 
               numeric(value = rep(curCoords[pointIter,dimIter], numWindows) >= lowerCoords[1:numWindows,dimIter] &
                          rep(curCoords[pointIter,dimIter], numWindows) < upperCoords[1:numWindows,dimIter], length = numWindows)
         }# dimIter
         density <- density + isInWindow
      }
   }# pointIter
   return(density)
})