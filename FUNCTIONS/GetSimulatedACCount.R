GetSimulatedACCount <- function( AC.sp, habitat.r)
   {
   r <- raster(extent(habitat.r))
   res(r) <- res(habitat.r)
   rtot <- rasterize(AC.sp@coords, r, fun="count")
   return(list("POP.r" = rtot))
   }
