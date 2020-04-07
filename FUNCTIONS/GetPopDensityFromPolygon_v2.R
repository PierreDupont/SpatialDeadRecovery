# GetPopDensityFromPolygon_v2 <- function(habPolygon = NULL,
#                               habResolution = NULL,
#                               habRaster = NULL,
#                               posterior.sigma = NULL,
#                               posterior.sxy,
#                               posterior.z,
#                               alive.states,
#                               plot.check){
#    
#    # Retrieve the resolution of the raster
#    if(!is.null(habRaster)){
#       habResolution <- res(habRaster)
#       }
#    
#    # Create a raster of the zone of interest
#    if(is.null(habRaster)){
#       habRaster <- raster(extent(habPolygon))
#       res(habRaster) <- habResolution
#       habRaster <- rasterize(habPolygon, habRaster)
#       }
#    
#    i<-1
#    id<-1
#    
#    
#    sxy.<-lapply(1:dim(posterior.sxy)[1],function(i)posterior.sxy[i,,])
#    sxy.<-do.call(rbind,sxy.)
#    dim(sxy.)
#    
#    z.<-lapply(1:dim(posterior.z)[1],function(i)posterior.z[i,])
#    z.<-do.call(c,z.)
#    dim(z.)
#    
#    sigma.<-rep(posterior.sigma,each=dim(posterior.z)[1]*dim(posterior.sxy)[2])
#    length(sigma.)
# 
#    keep<-z.%in%alive.states
#    sigma.<-sigma.[keep] 
#    z.<-z.[keep] 
#    sxy.<-sxy.[keep,] 
#    
#    
#    ACxy<-data.frame(x=sxy.[,1],y=sxy.[,2])
#    #ACxy$value<-1
#    coordinates(ACxy)<-ACxy
#    proj4string(ACxy)<-proj4string(habRaster)
#    ACxy@data
#    ACxy.sf <- st_as_sf(ACxy)
#    
#    habxy <- rasterToPoints(habRaster, spatial=TRUE)
#    habxy.sf <- st_as_sf(habxy)
#    
# 
#    
#    r<-raster::rasterize(coordinates(ACxy),habRaster)
#    r<-r/sum(r[],na.rm=TRUE)
#    r<-r*length(z.)/dim(posterior.z)[1]
#    sum(r[],na.rm=TRUE)
#    plot(r)
#    
#    
#    # temp<-st_buffer(ACxy.sf,0.1)
#    # attributes(temp)
#    #    temp<-fasterize(temp,habRaster,field = "value",fun="max")
#    # 
#    #    plot(temp)
# i<-1
# sigma2.2<-2*(sigma.*1000) ^2
#    out<-lapply(1:length(habRaster),function(i){
#       print(i)
#       dist <- as.numeric(st_distance(habxy.sf[i,],  ACxy.sf))
#       p0<- exp(-(dist^2)/sigma2.2)
#       sum(p0)
#    })
#    
#    dist<-st_distance(habxy.sf,  ACxy.sf)
#    
#   rr.sum.list<- lapply(1:dim(posterior.sxy)[1],function(i){
#      print(i)
#    #dim(posterio|r.sxy)
#    AC.xy<-posterior.sxy[i,posterior.z[i,]%in%alive.states,]
#    r.list<-   lapply(1:dim(AC.xy)[1],function(id){
#    r <- distanceFromPoints(habRaster,  posterior.sxy[i,id,])
#    #range(r[])
#    #plot(AC.xy[,])
#    #plot(r,add=TRUE)
#    #points(AC.xy[id,2]~AC.xy[id,1],col="red",pch=19,cex=1)
# 
#    
#    #r<-calc(r,function(x)exp(-(r[]^2)/(2*(posterior.sigma[i]*1000)^2)))
#    r[]<-exp(-(r[]^2)/(2*(posterior.sigma[i]*1000)^2))
#    r[]<-r[]/sum(r[])
#    #r[]<-dnorm( r[],sd=(posterior.sigma[i]*1000)^2)
#    #print(range(r[]))
#    return(r)
# })
#    r<-brick (r.list)
#    r.sum<-calc(r,sum)
#    return(r.sum)
#    })
#    
# #plot(r.list[[2]]) 
# 
# #plot(r.sum)
# #sum(r.sum[]) == dim(AC.xy)[1]
# 
#    #plot(distance.r)
#    #points(posterior.sxy[i,,],pch=19,cex=0.2)
#    #distance.r
# 
#    
#    # r <- raster(ncol=36,nrow=18)
#    # xy <- rbind(c(0,0),c(50,50))
#    # d1 <- distanceFromPoints(r, xy) 
#    # plot(d1)
#    # points(xy,pch=19,cex=0.2)
#    
#    
#    # JOE'S  FANTASTIC BULLSHIT
#    # Filter out posterior sxy for dead individuals
#    deadId <- which(posterior.z != alive.states, arr.ind = TRUE)
#    # for(i in 1:dim(deadId)[1]){
#    #    posterior.sxy[deadId[i,1],deadId[i,1],1:2] <-
#    # }
#    post.x <- posterior.sxy[ , ,1]
#    post.x[deadId] <- -999
#    posterior.sxy[,,1] <- post.x
#    
#    # Convert habRaster to SpatialGrid object 
#    spRaster <- as(habRaster, "SpatialGridDataFrame")
#    
#    # Calculate cell-specific density for each iteration
#    curCoords<-posterior.sxy[1,,]
#    Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
#       # Get the cell index for each coordinate
#       curGridIndeces <- getGridIndex(curCoords, getGridTopology(inRaster), all.inside = FALSE)
#       # Get rid of coordinate indeces that do not fall in a cell
#       curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
#       # Create a frequency table of grid indeces
#       freqTable <- table(as.character(curGridIndeces))
#       # Initialise an output vector
#       outVector <- rep(0, nrow(coordinates(inRaster)))
#       outVector[as.integer(names(freqTable))] <- freqTable
#       outVector[is.na(spRaster@data$layer)] <- NA
#       outVector
#    }, inRaster = spRaster, MARGIN = 1)
# 
#    # EXTRACT SUMMARY STATISTICS PER HABITAT CELL
#    meanDensity <- apply(Density, 1, mean)
#    CI <- apply(Density, 1, function(x)quantile(x, c(0.025,0.5,0.975), na.rm = TRUE))
#    
#    # FEED IN RASTERS
#    meanDensity.r <- lowCIDensity.r <- upCIDensity.r <- medianDensity.r <- habRaster
#    meanDensity.r[] <- meanDensity
#    lowCIDensity.r[] <- CI[1, ]
#    medianDensity.r[] <- CI[2, ]
#    upCIDensity.r[] <- CI[3, ]
#    
#    if(plot.check){
#       par(mfrow = c(1,2))
#       ## Plot mean density 
#       max.cut <- max(meanDensity.r[], na.rm = TRUE)
#       cuts <- seq(0, max.cut, length.out = 101) 
#       pal <- rev(terrain.colors(length(cuts)))
#       plot(meanDensity.r, breaks = cuts, col = pal, main ="Mean density", legend=FALSE, axes=FALSE)
#       plot(meanDensity.r, legend.only=TRUE, col=pal,
#            legend.width=1, legend.shrink=1,
#            axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
#                           labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
#                           cex.axis=0.6),
#            legend.args = list(text=paste("Density (ind.",(habResolution/1000)^2,"km-2)",sep=""), side=4, font=2, line=2.5, cex=0.8))
#       plot(habPolygon, add = TRUE)
#       
#       ## Plot CI width
#       max.cut <- max(upCIDensity.r[]-lowCIDensity.r[], na.rm = TRUE)
#       cuts <- seq(0, max.cut, length.out = 101) 
#       pal <- rev(heat.colors(length(cuts)))
#       plot(upCIDensity.r-lowCIDensity.r, breaks = cuts, col = pal, main ="95% CI width", legend = FALSE, axes=FALSE)
#       plot(upCIDensity.r, legend.only=TRUE, col=pal,
#            legend.width=1, legend.shrink=1,
#            axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
#                           labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
#                           cex.axis=0.6))
#       plot(habPolygon, add = TRUE)
#    }
#    
#    
#    return(list("mean.Density" = meanDensity.r,
#                "median.Density" = medianDensity.r,
#                "upperCI.Density" = upCIDensity.r,
#                "lowerCI.Density" = lowCIDensity.r))
# }
