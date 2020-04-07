#' @title Function to simulate detection of simulated individuals within an SCR framework, under the influence of spatial covariates (at the detector level).
#'
#' @description
#' \code{SimulateDetectionPoisson} simulates detection from a Poisson distribution. Then it can subsample detections from the Poisson and 
#'    do some pooling over larger grid size. it can pool for 3 types of model (Poisson, Bernoulli, Binomial). It returns a list of object similar to what the "normal" \code{SimulateDetection} fonction can do 
#' 
#' @param lambda0 Numeric variable denoting the intercept value of the negative binomial detection function describing the decay of detection with increasing distance from the AC.
#' @param sigma Numeric variable denoting the scale parameter of the detection function.
#' @param AC.sp Spatial points dataframe with individual activity center locations.
#' @param detector.sp Spatial points dataframe with detector locations. this should come from from MakeSearchGrid
#' @param detector.covs A data frame with the covariate values associated with each habitat grid cell (order identical to "covs").
#' @param betas A numeric vector denoting the coefficient values associated with the covariates to be used during simulations (effect on detection probability).
#' @param n.samples A numeric variable denoting the total number of samples to be returned.
#' @param alive.ids A numeric vector denoting the ids of the individuals that are alive. they should match the AC.sp. 
#' @param alpha A numeric variable denoting the probability that a sample is designated to an individual (i.e. implementation of partial identification)
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' @param seed A seed number for result to be reproducible (default is NULL)
#' areas on.
#' 
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SimulateDetection.R
#' @keywords simul
#'
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#'   det.sim <- SimulateDetectionPoisson(  lambda0 = 0.8
#'                                        , sigma = 6000
#'                                        , AC.sp = mySimulatedACs
#'                                        , detector.sp =  mySubDetectors$detector.sp
#'                                        , detector.covs = NULL
#'                                        , betas=NULL
#'                                        , seed = 10
#'                                        , n.samples = NULL
#'                                        , alive.ids = c(1:50) 
#'                                        , alpha = NULL
#'                                        , type = "Poisson"
#'                                        , plot.check = TRUE)






SimulateDetectionPoisson <- function( lambda0
                                     , sigma
                                     , AC.sp
                                     , detector.sp
                                     , detector.covs = NULL
                                     , betas
                                     , n.samples = NULL
                                     , alive.ids = 1:length(AC.sp)
                                     , alpha = NULL
                                     , type = "Poisson"
                                     , plot.check = TRUE
                                     , seed=  NULL){


   

      projection(AC.sp) <- projection(detector.sp)#---Just in case...
 
      ## log transform for covariates ## 
      fixed.effects <- rep(log(lambda0), length(detector.sp))
   
   
      #---- If covariates ----
   if(!is.null(detector.covs)){
      temp <- formula(paste("~", paste(names(detector.covs), collapse="+"), sep=""))
      Xmat <- model.matrix(temp, detector.covs)
      
      fixed.effects <- Xmat[,] %*% betas
      
      fixed.effects <- do.call(rbind,lapply(1:length(AC.sp), function(x) t(fixed.effects)))#---because the environmental covariates only change for different detectors, not for different ACs
      X <- data.frame(Xmat)
   }
   
      #---- Calculate distances between ACS and detectors ----
      #occasionally rownames pose a problem; fix now:
      dimnames(detector.sp@data)[[1]] <- dimnames(detector.sp@coords)[[1]] <- 1:length(detector.sp)
      dimnames(AC.sp@data)[[1]] <- dimnames(AC.sp@coords)[[1]] <- 1:length(AC.sp)

      D <- gDistance(detector.sp, AC.sp, byid=TRUE)
   
      #----  Compute the half normal    
      lambda0 <- exp(fixed.effects)
      lambda <- lambda0 * exp(-D*D/(2*sigma*sigma))
      
     
      #----  Sample observations  using a poisson process ----   
      #----  Set seed ----
      if(!is.null(seed)){
         y.seed <- D
         set.seed(seed)
         y.seed[] <- sample(10000, dim(D)[1]*dim(D)[2], replace=TRUE)
         temp <- abind(lambda, y.seed, along=3)
         y <- t(apply(temp, c(2,1), function(x){
              set.seed(x[2])
               rpois(1, x[1])}))
      }else{
         y <- t(apply(lambda, c(2,1), function(x) rpois(1, x)))}
         y.origin <- y
      
      # #if there are no ids alive, return a 0 matrix
      if(is.null(alive.ids) ){ alive.ids <- c(1:length(AC.sp))} 
      if(alive.ids[1]== 0 ){
         y[] <- 0    
      }else{
         y <- y[alive.ids, ]
      }
         
      #---- If subsampling ----
       
      if(!is.null(n.samples)){
         y.sub <- y     
         # if individual not alive then dont do the sampling on them
         if(length(alive.ids) != length(AC.sp)){ 
            y[ which(!(as.numeric(row.names(y)) %in% alive.ids)), ] <- 0} ### so the number of samples are not drawn from individuals that are not alive. 
         
         if( n.samples < sum(y.sub)){
            y.sub[] <-   apply(rmultinom(n.samples, 1, p=y/sum(y)), 1, sum)
            lambda0 <- lambda0 * sum(y.sub)/sum(y)#---to output correct parameter after subsampling
            betas[1] <- log(lambda0)[1] #---to output correct parameter after subsampling
            y <- y.sub
         }
        }
  
      
      y.orig <- y
      
      

      
      
      ## if there's no pooling then, have the possibility to transform it as a Bernoulli ## 
      if( sum(detector.sp$Id %in% detector.sp$main.cell.id[1:50])==50  ){ # do it over the 50 detect to save time ##
         if(type=="Bernoulli"){
            y[y>0] <- 1
         }
         
      }else{ # IF POOLING, pool Poisson, pool Binomial, pool Bernoulli ##
          

         
         #----  IF BINOMIAL ----
         if(type=="Binomial"){
            
            # Put everything to bernoulli ##    
            y.binom <- y
            y.binom[y>0] <- 1
           
            #ugly trick so when number of alive id== 1 it still works...
            if(dim(as.matrix(y))[2]==1 ){
            y.binom <- t(as.matrix(y.binom))
            } 
            
            # Aggregate over main detectors #
            y.binom <- t(aggregate(t(y.binom), by=list(detector.sp@data$main.cell.id), FUN=sum))
            row.names(y.binom) <- 1:dim(y.binom)[1]
            y.binom <- y.binom[-c(1), ]
            rownames(y.binom) <- 1:dim(y.binom)[1]
            y <- y.binom
            
         }else{## if Poisson and Bernoulli ## 
            if(type=="Poisson" | type=="Bernoulli"){
               # Aggregate over main detectors #
               #ugly trick so when number of alive id== 1 it still works...
               if(dim(as.matrix(y))[2]==1 ){
                  y <- t(as.matrix(y))
               } 
               
               y.pois <- t(aggregate(t(y), by=list(detector.sp@data$main.cell.id), FUN=sum))
               row.names(y.pois) <- 1:dim(y.pois)[1]
               y.pois <- y.pois[-c(1), ]
               rownames(y.pois) <- 1:dim(y.pois)[1]
               y <- y.pois
               ## if Bernoulli, convert from Poisson to Bern ##
               if(type=="Bernoulli"){ y[y>0] <- 1}

               
            }
         }
      }
      
      if(is.null(dim(y))) {y <- t(as.matrix(y))}
      
      #---- get number of trials ----
      y.trials <- unlist(lapply(sort(unique(detector.sp@data$main.cell.id)), function(id){
         sum(detector.sp@data$main.cell.id==id)
      }))
     
      
      
      sub.detector.sp <- detector.sp
      temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
      temp <- temp[order(temp$main.cell.id), ]## order by IDS
      
      main.detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], 
                                            data=data.frame(temp), proj4string=CRS(projection(detector.sp)))
      
      #------------------------------------------#
      #---- Plot Checks ----
      #------------------------------------------#
             AC.sp.alive <- AC.sp[alive.ids, ]
      
            if(plot.check == TRUE){
               par(mfrow=c(1,3))
     
               if(!is.null(dim(y.orig))) {
               detect <- which(rowSums(y.orig) > 0)
               set.seed(Sys.time())  ## cancel the set.seed 
               id <- sample(detect ,1)
               }else{y.orig <- t(as.matrix(y.orig))
                     id <- 1}
               
               #plot 1
               plot(sub.detector.sp, cex=0.1 , pch=16, main="Original Poisson", col="grey")
               points(sub.detector.sp, cex=DoScale(y.orig[id,],l = 0, 1.5), pch=16)
               points(AC.sp.alive[id,], pch=16, col="red")
               
               
               if(sum(detector.sp$Id %in% detector.sp$main.cell.id[1:50])!=50 ){
                  
               #plot 2
               plot(sub.detector.sp, cex=0.1 , pch=16, main="After pooling", col="grey")
               points(main.detector.sp, pch=16, col="black", cex=DoScale(y[id,],l = 0, 1.5))
               points(AC.sp.alive[id,], pch=16, col="red")
               
               #plot 3
               plot(main.detector.sp, col=grey(0.4), cex=0.6 , pch=16, main="Individual detections")

               #only keep alive ids 
               for( i in 1:length(AC.sp.alive)){#lapply( 1: length(AC.sp.alive), function(x) {
                  this.row <- y[i,]
                  this.det <- main.detector.sp[this.row>0, ]
                  if(length(this.det) > 0){
                     segments(coordinates(this.det)[ ,1], coordinates(this.det)[ ,2], coordinates(AC.sp.alive[i, ])[ ,1],
                              coordinates(AC.sp.alive[i, ])[ ,2], col="pink")
                  }
               }
               points(AC.sp.alive, col="red", pch=19)
               
               
               
               }else{
                  #plot 2
                  plot(sub.detector.sp, pch=16, col="grey", cex=DoScale(y[id,],l = 0.4, 1.5), main="Final Y")
                  points(AC.sp.alive[id,], pch=16, col="red")
                  
                  #plot 3
                  plot(sub.detector.sp, col="blue", pch=16, cex=0.3, main="Individual detections")
                  plot(sub.detector.sp, col="blue", pch=16, cex=0.3, add=TRUE)
                  
                  #only keep alive ids 
                  lapply( 1: length(AC.sp.alive), function(x) {
                     this.row <- y[x,]
                     this.det <- sub.detector.sp[this.row>0, ]
                     if(length(this.det) > 0){
                        segments(coordinates(this.det)[ ,1], coordinates(this.det)[ ,2], coordinates(AC.sp.alive[x, ])[ ,1],
                                 coordinates(AC.sp.alive[x, ])[ ,2], col="pink")
                     }
                  })
                  points(AC.sp.alive, col="red", pch=19)
                  
                  } 
               }
               

      
      #---- SLIM TO INDIVIDUALS DETECTED ---- 
      y.all <- matrix(0, nrow=length(AC.sp), ncol=dim(y)[2])
      y.all[alive.ids,] <- y
      detected <- apply(y, 1, max)>0
      y <- y[detected, ]
             
      
      #---- CONVERT DECTECTION HISTORY TO DETECTION EVENTS ---- 
      # depending if there's pooling or not, choose the good detectors 
      if(sum(detector.sp$Id %in% detector.sp$main.cell.id[1:50])!=50  ){
         detector <- main.detector.sp
      }else{
         detector <- sub.detector.sp
      }
      
      if(sum(y)==0   ){## if there are no detections, do return an empty det.sp 
         det.sp <- as(detector, "SpatialPoints")
         det.sp <- det.sp[-c(1:length(det.sp)),] 
      }else{
         if(is.null(dim(y))) {
            y <- t(as.matrix(y))
         }
         
         dimnames(y) <- list(1:dim(y)[1], 1:length(detector))  #---DO NOT REMOVE! :)
         
 
         
         raw.data <- do.call(rbind,lapply(1:dim(y)[1], function(x){
            
            this.row <- y[x,]
            these.det <- coordinates(detector)[rep(1:length(detector), this.row), ]
            if(is.null(dim(these.det))) these.det <- t(these.det)
            if(length(these.det)>0) data.frame(id=x,these.det)
            if(!is.null(dim(these.det))){ if(dim(these.det)<1) these.det <- data.frame(NA,NA)}
            dimnames(these.det)[[2]] <- c("x","y")
            return(na.omit(these.det))
         }))
         raw.data <- data.frame( raw.data)
         rownames(raw.data) <- 1:dim(raw.data)[1]
         det.sp <- SpatialPointsDataFrame(raw.data[,c("x","y")], data=raw.data, proj4string=CRS(projection(detector.sp)))#
         
         dimnames(AC.sp@coords) <- list(c(1:length(AC.sp)), c("x","y"))
         if(plot.check==T){points(det.sp,cex=0.8)}
      }

        
      #----DO THE PARTIAL ID IMPLEMENTATION ----
      if(type=="Poisson" & !is.null(alpha)){
         ind.y <- y
         
         ind.y[] <- unlist(lapply(y, function(x) rbinom(1, x, alpha)))
         ind.y <- ind.y[apply(ind.y, 1, function(x) max(x)>0), ]#slim it down to detected individuals
         count.y <- apply(y,2,sum) #---everything: known and unknown ids
         
         y <- ind.y
      }
      
      
      
      out<-list(
           y.origin=y.origin #This is the original y (poisson) before aggregation. 
         , y = y         #---this is the y for analysis
         , y.all = y.all #---this is the y to match with simulated ACs (contains all individuals, including those that were never detected)
         , D = D
         , lambda0 = lambda0
         , sigma = sigma
         , betas = betas
         , detector.covs = detector.covs
         , det.sp = det.sp
         , alpha = alpha # probability of making an ID
         , n.trials = y.trials
      )  
      
      if(type=="Poisson" & !is.null(alpha))out$count.y<-count.y
      
      return(out)
               
}
            
            
            

      

   
   