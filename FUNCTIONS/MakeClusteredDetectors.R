MakeClusteredDetectors <- function(  raster.origin
                                   ,  detector.sp
                                   , fact= 2
                                   , even.col= TRUE
                                   , plot.check=TRUE){
   
   raster.origin1 <- raster.origin
   
   if(fact>1){
    r <- aggregate(raster.origin1, fact=fact, expand=TRUE)
   }else{
      r <- raster.origin
   }
   
   ##
    r[] <- 1:ncell(r)
    id.det.cell <- cellFromXY(r, detector.sp)
    
    ###
   
   temp <- as.matrix(r)
   out <- list()
   
   
   odd.col <- (1:dim(temp)[2]) %% 2
   if(even.col==TRUE){
      odd.col <- as.numeric(odd.col == 1)
   }else{
      odd.col <- as.numeric(odd.col !=1)
   }
   

   for(i in 1:dim(temp)[1]){
      out[[i]] <- temp[i, odd.col==1 -(i %% 2) ] 
   }
   
   out <- do.call(c, out)
   out <- out[!is.na(out)]
   
   if(plot.check==TRUE){
   plot(raster.origin)
   points(detector.sp)
   points(detector.sp[which(id.det.cell %in% out),], col="red", pch=16)
   }
   
   ID <- which(id.det.cell %in% out)
   
   detector.sp <- detector.sp[which(id.det.cell %in% out),]
   
   return(list(  main.cell.id = detector.sp$main.cell.id
               ,Id = detector.sp$Id
               ,ID.pos.object = ID))
}