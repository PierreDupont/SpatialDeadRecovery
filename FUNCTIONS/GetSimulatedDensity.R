GetSimulatedDensity <- function( AC.sp = mySimulatedACs
                            , habitat.r = myHabitat$habitat.r
                            , sigma = 1)
   {
   ##  Calculate distances between each habitat cells and individual ACs 
   habitat.xy <- coordinates(habitat.r)
   sx <- AC.sp$x
   sy <- AC.sp$y
   D2 <- apply(habitat.xy, 1, function(c){(sy-c[2])^2 + (sx-c[1])^2})

   ## Calculate INDIVIDUAL SPACE-USE (i.e. percentage of time spent in a given habitat cell)
   lambda0 <- exp(-D2/(2*sigma^2))
   lambda <- lambda0/rowSums(lambda0)
   
   ## Convert INDIVIDUAL SPACE-USE matrices to a list of raster
   IndividualSU.list <- lapply(1:dim(lambda)[1], function(x){
      r <- habitat.r
      raster::values(r) <- lambda[x, ]
      return(r)
      })
      
   ## Store INDIVIDUAL SPACE-USE rasters in a raster-stack
   IndividualSU.r.stack <- stack(IndividualSU.list)

   ## Sum the raster stack to get POPULATION SPACE-USE 
   return(list( "ID.r" = IndividualSU.r.stack
              , "POP.r" = sum(IndividualSU.r.stack)))
   }


