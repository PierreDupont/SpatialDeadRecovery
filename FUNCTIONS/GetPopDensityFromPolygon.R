GetPopDensityFromPolygon <- function(habPolygon = NULL,
                              habResolution = NULL,
                              habRaster = NULL,
                              posterior.sxy,
                              posterior.z,
                              alive.states,
                              plot.check){
   
   # Retrieve the resolution of the raster
   if(!is.null(habRaster)){
      habResolution <- res(habRaster)
      }
   
   # Create a raster of the zone of interest
   if(is.null(habRaster)){
      habRaster <- raster(extent(habPolygon))
      res(habRaster) <- habResolution
      habRaster <- rasterize(habPolygon, habRaster)
      }
   
   # JOE'S  FANTASTIC BULLSHIT
   # Filter out posterior sxy for dead individuals
#<<<<<<< HEAD

   #---[RB] fixed 2019-08-29
   deadId <- which(!posterior.z %in% alive.states)#, arr.ind = TRUE)
   posterior.sxy[,,1][deadId]<--999
   #hist(apply(posterior.sxy[,,1],1,function(x)sum(x!=-999)))
   
  
# =======
#    deadId <- which(posterior.z != alive.states, arr.ind = TRUE)
#    post.x <- posterior.sxy[ , ,1]
#    post.x[deadId] <- -999
#    posterior.sxy[,,1] <- post.x
# >>>>>>> 8de688948e7224b9112cb5768978cb5f5d0646d7
   
   # Convert habRaster to SpatialGrid object 
   spRaster <- as(habRaster, "SpatialGridDataFrame")
   
   # Calculate cell-specific density for each iteration
   Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
      # Get the cell index for each coordinate
      curGridIndeces <- getGridIndex(curCoords, getGridTopology(inRaster), all.inside = FALSE)
      # Get rid of coordinate indeces that do not fall in a cell
      curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
      # Create a frequency table of grid indeces
      freqTable <- table(as.character(curGridIndeces))
      # Initialise an output vector
      outVector <- rep(0, nrow(coordinates(inRaster)))
      outVector[as.integer(names(freqTable))] <- freqTable
      outVector[is.na(spRaster@data$layer)] <- NA
      outVector
   }, inRaster = spRaster, MARGIN = 1)

   # EXTRACT SUMMARY STATISTICS PER HABITAT CELL
   meanDensity <- apply(Density, 1, mean)
   CI <- apply(Density, 1, function(x)quantile(x, c(0.025,0.5,0.975), na.rm = TRUE))
   
   # FEED IN RASTERS
   meanDensity.r <- lowCIDensity.r <- upCIDensity.r <- medianDensity.r <- habRaster
   meanDensity.r[] <- meanDensity
   lowCIDensity.r[] <- CI[1, ]
   medianDensity.r[] <- CI[2, ]
   upCIDensity.r[] <- CI[3, ]
   
   if(plot.check){
      par(mfrow = c(1,2))
      ## Plot mean density 
      max.cut <- max(meanDensity.r[], na.rm = TRUE)
      cuts <- seq(0, max.cut, length.out = 101) 
      pal <- rev(terrain.colors(length(cuts)))
      plot(meanDensity.r, breaks = cuts, col = pal, main ="Mean density", legend=FALSE, axes=FALSE)
      plot(meanDensity.r, legend.only=TRUE, col=pal,
           legend.width=1, legend.shrink=1,
           axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                          labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                          cex.axis=0.6),
           legend.args = list(text=paste("Density (ind.",(habResolution/1000)^2,"km-2)",sep=""), side=4, font=2, line=2.5, cex=0.8))
      plot(habPolygon, add = TRUE)
      
      ## Plot CI width
      max.cut <- max(upCIDensity.r[]-lowCIDensity.r[], na.rm = TRUE)
      cuts <- seq(0, max.cut, length.out = 101) 
      pal <- rev(heat.colors(length(cuts)))
      plot(upCIDensity.r-lowCIDensity.r, breaks = cuts, col = pal, main ="95% CI width", legend = FALSE, axes=FALSE)
      plot(upCIDensity.r, legend.only=TRUE, col=pal,
           legend.width=1, legend.shrink=1,
           axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                          labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                          cex.axis=0.6))
      plot(habPolygon, add = TRUE)
   }
   return(list("mean.Density" = meanDensity.r,
               "median.Density" = medianDensity.r,
               "upperCI.Density" = upCIDensity.r,
               "lowerCI.Density" = lowCIDensity.r))
}
