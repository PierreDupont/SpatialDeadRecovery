GetPredictedDensity <- function( posterior.sxy = my.jags.output$sims.list$sxy
                               , posterior.z = my.jags.output$sims.list$z
                               , posterior.sigma = my.jags.output$sims.list$sigma
                               , habitat.r = myHabitat$habitat.r
                               , habitat.mx = myHabitat$habitat.mx
                               , alive.states = NULL)
{
   ## SET VARIABLES
   n.iter <- dim(posterior.z)[1]
   n.individuals <- dim(posterior.z)[2]
   n.years <- ifelse(length(dim(posterior.z)) == 2, 1, dim(posterior.z)[3])
   n.cells <- length(which(raster::values(habitat.r)==1))

   if(is.null(alive.states)){alive.states <- 1}   
   
   ## REFORMAT INPUT OBJECTS
   sx <- sy <- Z <- z <- sigma <- array(0, dim = c(n.iter, n.individuals, n.years)) 
   if(length(dim(posterior.sxy)) == 4){
      sy <- posterior.sxy[ , ,2, ]
      sx <- posterior.sxy[ , ,1, ]
      }#if
   if(length(dim(posterior.sxy)) == 3){
      sy[ , ,1:n.years] <- posterior.sxy[ , ,2]
      sx[ , ,1:n.years] <- posterior.sxy[ , ,1]
      }#if
   z[] <- posterior.z
   sigma[] <- posterior.sigma
   
   ## Subset habitat.xy to match active cells.
   habitat.xy <- which(t(habitat.mx) == 1, arr.ind = TRUE)

   ## Keep AC locations, Iterations and time of alive individuals only (z == 1)
   Z[z %in% alive.states] <- 1
   SY <- sy[Z==1]
   SX <- sx[Z==1]
   SIGMA <- sigma[Z==1]
   TIME <- which(Z == 1, arr.ind = TRUE)[ ,3] 
   ITER <- which(Z == 1, arr.ind = TRUE)[ ,1] 
   
   ## Calculate distance between ACs and habitat cells
   print("Calculating individual AC-Habitat DISTANCE MATRIX :")
   pb <- txtProgressBar(min = 0, max = length(habitat.xy[,1]), style = 3)
   D2 <- apply(habitat.xy, 1, function(c){
      temp.D2 <- (SY-c[2])^2 + (SX-c[1])^2
      setTxtProgressBar(pb, which(habitat.xy[,1] == as.numeric(c[1]) & habitat.xy[,2] == as.numeric(c[2])))
      return(temp.D2)
      })
   close(pb)
   
   ## Create arrays of Population counts and density 
   SPACE.USE <- array(0, c(n.cells, n.iter, n.years))
   print("Extracting individual SPACE-USE RASTERS :")
   pb <- txtProgressBar(min = 0, max = length(SX), style = 3)
   for(x in 1:length(SX)){
      curMap <- exp(-D2[x, ]/(2*SIGMA[x]^2))
      curMap <- curMap/(sum(curMap))                                   
      SPACE.USE[ , ITER[x], TIME[x]] <- SPACE.USE[ , ITER[x], TIME[x]] + curMap
      setTxtProgressBar(pb, x)
      }#x
   close(pb)

   ## Create raster stacks of yearly population densities and individuals counts
   for(t in 1:n.years){
      mean.SU.r <- q2.5.SU.r <- q97.5.SU.r <- habitat.r
      mean.SU.r[] <- q2.5.SU.r[] <- q97.5.SU.r[] <- 0
      
      raster::values(mean.SU.r)[which(habitat.r[] == 1)] <- apply(SPACE.USE[ , ,t], 1, mean)
      raster::values(q2.5.SU.r)[which(habitat.r[] == 1)] <- apply(SPACE.USE[ , ,t], 1, function(x)quantile(x,.025))
      raster::values(q97.5.SU.r)[which(habitat.r[] == 1)] <- apply(SPACE.USE[ , ,t], 1, function(x)quantile(x,.975))
      
      if(t>1){
         mean.rstack <- stack(mean.rstack, mean.SU.r)
         q2.5.rstack <- stack(q2.5.rstack, q2.5.SU.r)
         q97.5.rstack <- stack(q97.5.rstack, q97.5.SU.r)
         }else{
         mean.rstack <- mean.SU.r
         q2.5.rstack <- q2.5.SU.r
         q97.5.rstack <- q97.5.SU.r
         }#if
      }#t
   
   ## Create Output list
   return(list( "mean" = mean.rstack
              , "q2.5" = q2.5.rstack
              , "q97.5" = q97.5.rstack))
}

