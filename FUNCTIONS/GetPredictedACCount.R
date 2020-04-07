GetPredictedACCount <- function( posterior.sxy = my.jags.output$sims.list$sxy
                               , posterior.z = my.jags.output$sims.list$z
                               , habitat.r = myHabitat$habitat.r
                               , habitat.mx = myHabitat$habitat.mx
                               , alive.states = NULL)
{
   ## SET VARIABLES
   n.iter <- dim(posterior.z)[1]
   n.individuals <- dim(posterior.z)[2]
   n.years <- ifelse(length(dim(posterior.z)) == 2, 1, dim(posterior.z)[3])
   n.cells <- sum(raster::values(habitat.r) == 1)
   if(is.null(alive.states)){alive.states <- 1}  
   
   ## REFORMAT INPUT OBJECTS
   sx <- sy <- Z <- z <- sigma <- array(0, dim = c(n.iter, n.individuals, n.years)) 
   if(length(dim(posterior.sxy)) == 4){
      sy <- posterior.sxy[ , ,2, ]
      sx <- posterior.sxy[ , ,1, ]
      }#if
   if(length(dim(posterior.sxy)) == 3){
      sy[,,1:n.years] <- posterior.sxy[ , ,2]
      sx[,,1:n.years] <- posterior.sxy[ , ,1]
      }#if
   z[] <- posterior.z
   y.max <- dim(habitat.mx)[1]
   x.max <- dim(habitat.mx)[2]
   habitat.xy <- expand.grid(1:x.max,1:y.max)
   
   ## Keep AC locations, Iterations and time of alive individuals only (z == 1)
   Z[which(z %in% alive.states)] <- 1
   SY <- sy[Z==1]
   SX <- sx[Z==1]
   SX <- trunc(SX)+1
   SY <- trunc(SY)+1
   TIME <- which(Z == 1, arr.ind = TRUE)[ ,3] 
   ITER <- which(Z == 1, arr.ind = TRUE)[ ,1] 

   ## Create arrays of Population counts and density 
   AC.Pop.Density <- array(0, c(y.max, x.max, n.iter, n.years))
   print("Extracting individual AC LOCATIONS :")
   pb <- txtProgressBar(min = 0, max = length(SX), style = 3)
   for(x in 1:length(SX)){
      AC.Pop.Density[SY[x],SX[x],ITER[x],TIME[x]] <- AC.Pop.Density[SY[x],SX[x],ITER[x],TIME[x]] + 1
      setTxtProgressBar(pb, x)
      }#x
   close(pb)
   
   ## Create raster stacks of yearly population densities and individuals counts
   for(t in 1:n.years){
      mean.DENS.r <- raster(apply(AC.Pop.Density[ , , ,t],c(1,2),mean))
      q2.5.DENS.r <- raster(apply(AC.Pop.Density[ , , ,t],c(1,2),function(x)quantile(x,.025)))
      q97.5.DENS.r <- raster(apply(AC.Pop.Density[ , , ,t],c(1,2),function(x)quantile(x,.975)))
      extent(mean.DENS.r) <- extent(q2.5.DENS.r) <- extent(q97.5.DENS.r) <- extent(habitat.r)
      res(mean.DENS.r) <- res(q2.5.DENS.r) <- res(q97.5.DENS.r) <- res(habitat.r)
      
      if(t>1){
         mean.rstack <- stack(mean.rstack, mean.DENS.r)
         q2.5.rstack <- stack(q2.5.rstack, q2.5.DENS.r)
         q97.5.rstack <- stack(q97.5.rstack, q97.5.DENS.r)
         }else{
         mean.rstack <- mean.DENS.r
         q2.5.rstack <- q2.5.DENS.r
         q97.5.rstack <- q97.5.DENS.r
         }#if
      }#t
   
   ## Create Output list
   return(list( "mean" = mean.rstack
              , "q2.5" = q2.5.rstack
              , "q97.5" = q97.5.rstack))
}