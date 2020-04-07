#' @title Function to simulate individual detections within an SCR framework, under the influence of 
#' different individual and spatial (detector-level) covariates.
#'
#' @description
#' \code{SimulateDetection} returns a list object with \code{} 
#' 
#' @param p0 \code{Numeric} variable denoting the intercept value of the half-normal detection function describing the 
#' decay of detection probability with increasing distance from the AC.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param AC.sp A \code{SpatialPointsDataFrame} object with individual activity center locations.
#' @param detector.sp  A \code{SpatialPointsDataFrame} object with detectors' locations.
#' @param id.covs.p0 A \code{DataFrame} with the individual-level covariates to affect p0.
#' @param id.betas.p0 A \code{vector} object with the coefficients associated to each individual covariate. 
#' @param det.covs.p0 A \code{DataFrame} with the detector-level covariates  to affect p0.
#' @param det.betas.p0 A \code{vector} object with the coefficients associated to each detector covariate. 
#' @param id.covs.sigma A \code{DataFrame} with the individual-level covariates to affect sigma. 
#' @param id.betas.sigma A \code{vector} object with the coefficients associated to each individual covariate. 
#' @param link.sigma A \code{character} to specify the link-function to be used when modelling the effect of covariates on sigma (default is exponential to avoid negative sigmas)
#' @param n.samples A \code{numeric} object denoting the total number of samples to be returned.
#' @param alive.ids A \code{vector} vector denoting the ids of the individuals that are alive. They should match the AC.sp. 
#' @param alpha A \code{probability} denoting fthe identification probability (implementation of partial identification for Poisson simulations)
#' @param type A \code{character} to specify the type of detections to be modelled; can be binary (\code{"Bernoulli"} or \code{"Binomial})
#' or it can be derived from a count process (\code{"Poisson"} or \code{"BinomFomPoisson}).
#' @param plot.check A \code{logical} for whether (\code{TRUE}) or not (\code{FALSE}) plots are to be generated during simulations.
#' @param seed A seed number for result to be reproducible (default is NULL)

#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' det.sim<-SimulateDetection(p0=0.6,sigma=1200,AC.sp=AC.sp,detector.sp=detector.sp,covs=covs, betas=betas,plot=TRUE)

SimulateDetectionAlpha <- function( p0 
                               , sigma
                               , AC.sp
                               , detector.sp
                               , id.covs.p0 = NULL
                               , id.betas.p0 = NULL
                               , det.covs.p0 = NULL
                               , det.betas.p0 = NULL
                               , id.covs.sigma = NULL
                               , id.betas.sigma = NULL
                               , link.sigma = "exponential"
                               , n.samples = NULL
                               , alive.ids = NULL
                               , alpha = NULL
                               , type = "Bernoulli"
                               , plot.check = TRUE
                               , habitat.poly = NULL
                               , seed = NULL
                               , typeD = 1)
   {
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== I.CLEAN & SET-UP THE INPUT DATA ====
   projection(AC.sp) <- projection(detector.sp)
   myP0 <- p0         ## For plotting purpose
   mySIGMA <- sigma   ## For plotting purpose
   
   #--- Create separate detectors and sub-detectors SpatialPointsDataFrame when needed
   if(type == "Binomial" | type == "Bernoulli")
      {
      #--- Subset detector.sp to main detectors only
      sub.detector.sp <- detector.sp
      temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
      temp <- temp[order(temp$main.cell.id), ]
      detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], data = data.frame(temp), proj4string = CRS(projection(detector.sp)))
      
      #--- Extract detector-specific number of trials
      n.trials <- unlist(lapply(detector.sp@data$main.cell.id,function(x){sum(sub.detector.sp@data$main.cell.id==x)}))
      n.trials <- matrix(n.trials, length(AC.sp), length(detector.sp), byrow = TRUE) 
      if(type == "Bernoulli"){n.trials[] <- 1}
      }
   
   #--- Check if detectors and covariates dimensions match
   if(!is.null(det.covs.p0)){
      if(dim(det.covs.p0)[1] != length(detector.sp)){
         stop("Detectors' covariates dataframe dimensions do not match detectors SpatialPointDataFrame dimensions!")                
         }
      }
   
   #--- Occasionally rownames pose a problem; fix now:
   dimnames(detector.sp@data)[[1]] <- dimnames(detector.sp@coords)[[1]] <- 1:length(detector.sp)
   dimnames(AC.sp@data)[[1]] <- dimnames(AC.sp@coords)[[1]] <- 1:length(AC.sp)
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== II.CALCULATE DETECTOR & INDIVIDUAL-SPECIFIC p0s ====
      ### ---- 1.Detector-specific covariates effects on p0 ----
   det.fixed.effects.p0 <- 0
   if(!is.null(det.covs.p0)){
      temp <- formula(paste("~", paste(names(det.covs.p0), collapse="+"),"-1", sep=""))
      Xmat <- model.matrix(temp, det.covs.p0)
      det.fixed.effects.p0 <- Xmat%*%det.betas.p0
      det.fixed.effects.p0 <- do.call(rbind,lapply(1:length(AC.sp),function(x)t(det.fixed.effects.p0)))
      X.det.p0 <- data.frame(Xmat)
      }
   
      ### ---- 2.Individual covariates effects on p0 ----
   id.fixed.effects.p0 <- 0
   if(!is.null(id.covs.p0)){
      temp <- formula(paste("~", paste(names(id.covs.p0),collapse="+"),"-1",sep=""))
      Xmat <- model.matrix(temp, id.covs.p0)
      id.fixed.effects.p0 <- Xmat%*%id.betas.p0
      id.fixed.effects.p0 <- do.call(cbind,lapply(1:length(detector.sp), function(x)id.fixed.effects.p0))
      X.ind.p0 <- data.frame(Xmat)
      }
   
      ### ---- 3.Calculate individual & detector-specific p0 ----
   if(type == "Bernoulli" | type == "Binomial"){
      intercept.p0 <- matrix(logit(p0),length(AC.sp),length(detector.sp))
      p0 <- inv.logit(intercept.p0 + id.fixed.effects.p0 + det.fixed.effects.p0)
      }
   if(type == "Poisson"){
      intercept.p0 <- matrix(log(p0), length(AC.sp), length(detector.sp))
      p0 <- exp(intercept.p0 + id.fixed.effects.p0 + det.fixed.effects.p0)
      }
   
   individual.p0<-p0
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== III.CALCULATE INDIVIDUAL-SPECIFIC sigmas ====
      ### ---- 1.Individual covariates effects on sigma ----
   id.fixed.effects.sigma <- 0
   if(!is.null(id.covs.sigma)){
      temp <- formula(paste("~", paste(names(id.covs.sigma),collapse="+"),"-1",sep=""))
      Xmat <- model.matrix(temp, id.covs.sigma)
      id.fixed.effects.sigma <- Xmat%*%id.betas.sigma
      X.ind.p0 <- data.frame(Xmat)
      }
      ### ---- 2.Calculate individual-specific sigma ----
   if(link.sigma=="exponential"){
      intercept.sigma <- rep(log(sigma),length(AC.sp))
      sigma <- exp(intercept.sigma + id.fixed.effects.sigma)
      }else{
      intercept.sigma <- rep(sigma,length(AC.sp))
      sigma <- intercept.sigma + id.fixed.effects.sigma
      }
   sigma <- matrix(sigma, length(AC.sp), length(detector.sp),byrow = FALSE)
   individual.sigma<-sigma
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== IV.CALCULATE INDIVIDUAL & DETECTOR-SPECIFIC p ====
   #--- Calculate distance matrix between detectors and ACs
   if(typeD == 1){D <- rgeos::gDistance(detector.sp, AC.sp, byid=TRUE)}
   if(typeD == 2){  
      x1 <- detector.sp$main.cell.x
      y1 <- detector.sp$main.cell.y
      x2 <- AC.sp$x
      y2 <- AC.sp$y
      D <- matrix(unlist(lapply(1:length(x2), function(dd){sqrt((x1 - x2[dd])^2 + (y1 - y2[dd])^2) })), length(x2), length(x1), byrow = TRUE)
      }#if  
 
   #--- Calculate Individual & Detector-specific detection probability
   P <- p0*exp(-D*D/(2*sigma*sigma))
   
   #--- Set dead individuals detection probability to 0
   if(!is.null(alive.ids)){
      P2 <- matrix(0, dim(P)[1], dim(P)[2])
      P2[alive.ids, ] <- P[alive.ids, ]
      P <- P2
      }
   
   ##-------------------------------------------------------------------------------------------------------------------
   ## ==== V.BERNOULLI or BINOMIAL DETECTION PROCESS ====
   if(type == "Binomial" | type == "Bernoulli")
      {      
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(n.trials), replace = TRUE)
      #-- Sample number of individual detections per main detector
      temp <- abind(n.trials, P, y.seed, along = 3)
      y <- apply(temp, c(1,2), function(x){
         set.seed(x[3])
         rbinom(1, x[1], x[2])})
      
      
      #--- Partial identification: sample individuals unidentified
      y.all.ind<-ind.y <- y
      if(!is.null(alpha)){
         ind.y[] <- unlist(lapply(y, function(x){rbinom(1, x, alpha)}))
         count.y <- apply(y, 2, sum)                                         # Counts data: known and unknown ids
         y <- ind.y                                                          # Identified IDs only
      }
      
      }              
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VI.POISSON DETECTION PROCESS ====
   if(type == "Poisson")
      {
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(P), replace = TRUE)
      
      #-- Sample number of samples found per detector
      temp <- abind(y.seed, P,along = 3)
      y <- apply(temp, c(1,2), function(x){
         set.seed(x[1])
         rpois(1, x[2])}) 
      
      #--- Resample to match n.samples per individual
      if(!is.null(n.samples)){
         if(n.samples < sum(y)){
            y.sub <- y       
            y.sub[] <- apply(rmultinom(n.samples, 1, p = y/sum(y)), 1, sum)
            p0 <- p0*sum(y.sub)/sum(y)                                        # to output correct parameter after subsampling
            y <- y.sub
            }
         }
      
      #--- Partial identification: sample individuals unidentified
      y.all.ind<-ind.y <- y
      if(!is.null(alpha)){
         ind.y[] <- unlist(lapply(y, function(x){rbinom(1, x, alpha)}))
         count.y <- apply(y, 2, sum)                                         # Counts data: known and unknown ids
         y <- ind.y                                                          # Identified IDs only
         }
      }
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VII.BINOMIAL FROM POISSON DETECTION PROCESS ====
   if(type == "BinomFromPoisson")
      { 
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(P), replace = TRUE)
         
      #-- Sample number of samples found per detector
      temp <- abind(y.seed, P,along = 3)
      y <- apply(temp, c(1,2), function(x){
               set.seed(x[1])
               rpois(1, x[2])}) 
         
      #--- Resample to match n.samples per individual
      if(!is.null(n.samples)){
         if(n.samples < sum(y)){
            y.sub <- y       
            y.sub[] <- apply(rmultinom(n.samples, 1, p = y/sum(y)), 1, sum)
            p0 <- p0*sum(y.sub)/sum(y)                                        # to output correct parameter after subsampling
            y <- y.sub
            }
         }
         
      #--- Partial identification: sample individuals unidentified
      y.all.ind<-ind.y <- y
 
      if(!is.null(alpha)){                              
         ind.y[] <- unlist(lapply(y, function(x){rbinom(1, x, alpha)}))
         count.y <- apply(y, 2, sum)                                         # Counts data: known and unknown ids
         y <- ind.y                                                          # Identified IDs only
         }
   
      #--- Aggregate Poisson counts to Bernoulli detections at the sub-detector level
      y.binom <- y
      y.binom[y>0] <- 1
      
      #--- Aggregate to main detectors
      y.binom <- t(aggregate(t(y.binom), by=list(detector.sp@data$main.cell.id), FUN=sum))
      row.names(y.binom) <- 1:dim(y.binom)[1]
      y.binom <- y.binom[-c(1), ]
      rownames(y.binom) <- 1:dim(y.binom)[1]
      y <- y.binom
      
      #--- Subset detector.sp to main detectors only
      sub.detector.sp <- detector.sp
      temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
      temp <- temp[order(temp$main.cell.id), ]
      detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], data=data.frame(temp), proj4string=CRS(projection(detector.sp)))
      }
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VIII.DEAD RECOVERY PROCESS ====
   if(type == "Recovery")
      { 
      added.P <- rep(1, dim(P)[1])
      added.P[alive.ids] <- 0
      P <- cbind(P,added.P)
      y <- apply(P, 1, function(x){which(rmultinom(1,1,x)==1)})
      }  
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VIII.SET-UP THE OUTPUT DATA ====
   # #--- MAYBE LATER...
   # raw.data<-do.call(rbind,lapply(1:dim(y)[1],function(x){
   #    
   #    this.row<-y[x,]
   #    these.det<-coordinates(detector.sp)[rep(1:length(detector.sp),this.row),]
   #    if(is.null(dim(these.det)))these.det<-t(these.det)
   #    if(length(these.det)>0) data.frame(id=x,these.det)
   #    if(!is.null(dim(these.det))){if(dim(these.det)<1)these.det<-data.frame(NA,NA)}
   #    dimnames(these.det)[[2]]<-c("x","y")
   #    return(na.omit(these.det))
   # }))
   # raw.data<-data.frame( raw.data)
   # rownames(raw.data)<-1:dim(raw.data)[1]
   # head(raw.data)
   # det.sp<-SpatialPointsDataFrame(raw.data[,c("x","y")],data=raw.data,proj4string=CRS(projection(detector.sp)))
   
   #--- Fix the names
   if(!is.vector(y)){dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])}  #---DO NOT REMOVE! :)
   dimnames(AC.sp@coords) <- list(c(1:length(AC.sp)), c("x","y"))
    
   #--- Slim to individuals detected
   y.all <- y
   if(length(dim(y)) > 1){
      detected <- apply(y, 1, max) > 0
      if(!any(detected)){detected <- rep(TRUE, dim(y)[1])}
      y <- y[detected, ]
      }
   if(length(dim(y)) == 0){
      detected <- which(y < dim(P)[2])
      if(!any(detected)){detected <- rep(TRUE, length(y))}
      y <- y[detected]
      }
   
   #--- Output list
   out <- list( y = y        
              , y.all = y.all
              , y.all.ind = y.all.ind
              , D = D
              , p0 = p0
              , sigma = sigma
              , id.covs.p0 = data.frame(id.covs.p0[detected, ])
              , id.betas.p0 = id.betas.p0
              , det.covs.p0 = det.covs.p0
              , det.betas.p0 = det.betas.p0
              , id.covs.sigma = data.frame(id.covs.sigma[detected, ])
              , id.betas.sigma = id.betas.sigma
              , n.samples = n.samples
              , alpha = alpha
              , individual.sigma=individual.sigma
              , individual.p0=individual.p0)  
   
   #if(type == "Poisson" & !is.null(alpha)){
      out$count.y <- count.y
   #   }
   if(type == "Binomial" | type == "BinomFromPoisson"){out$n.trials <- n.trials[1, ]}
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== IX.PLOT DETECTIONS ====
   if(plot.check)
      {
      ### ---- 1.Plot individual cov effects on p0 ----
      if(!is.null(id.covs.p0) & type != "BinomFromPoisson" & type != "Recovery")
         {
         AC.sp@data <- id.covs.p0
         for(i in 1:dim(id.covs.p0)[2])
            {
            myCov <- id.covs.p0[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(heat.colors(n = length(myValues), alpha = 1))
            
            par(mfrow=c(1,2))
            
            ## Plot the effect of covariate i on p0
            myBeta <- id.betas.p0[i]
            myX <- seq(min(myValues),max(myValues),length.out = 1000)
            if(type == "Poisson"){
               myY <- exp(log(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(id.covs.p0[i]), ylab = "lambda0", ylim = c(0,max(myY)), col = rev(heat.colors(length(myX))))
               }
            if(type != "Poisson"){
               myY <- inv.logit(logit(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(id.covs.p0[i]), ylab = "p0", ylim = c(0,1), col = rev(heat.colors(length(myX))))
               }
            
            ## Plot the Detections
            plot(AC.sp, main = names(myCov))
            plot(detector.sp, col="gray60", pch = 19, cex=0.6, add=TRUE)
            
            for(j in 1:length(myValues)){
               myIds <- which(AC.sp@data[ ,i] == myValues[j])
               plot(AC.sp[myIds, ], col = myCols[j], pch = 19, add = TRUE)
               
               lapply(myIds,function(x){
                  this.row <- y.all[x, ]
                  if(sum(this.row)>0){
                     this.det <- detector.sp[this.row>0, ]
                     segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                               , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                               , col = myCols[j])
                     }
                  })
               }
            }
         }
      
      ### ---- 2.Plot detector cov effects on p0 ----
      if(!is.null(det.covs.p0) & type != "BinomFromPoisson" & type != "Recovery")
      {
         detector.sp@data <- det.covs.p0
         for(i in 1:dim(det.covs.p0)[2])
         {
            myCov <- det.covs.p0[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(terrain.colors(n = length(myValues), alpha = 1))
            
            par(mfrow=c(1,2))
            
            ## Plot the effect of covariate i on p0
            myBeta <- det.betas.p0[i]
            myX <- seq(min(myValues),max(myValues),length.out = 1000)
            if(type == "Poisson"){
               myY <- exp(log(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(det.covs.p0[i]), ylab = "lambda0", ylim = c(0,max(myY)), col = rev(terrain.colors(length(myX))))
            }
            if(type != "Poisson"){
               myY <- inv.logit(logit(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(det.covs.p0[i]), ylab = "p0", ylim = c(0,1), col = rev(terrain.colors(length(myX))))
            }
            
            ## Plot the Detections
            plot(AC.sp)
            
            for(j in 1:length(myValues)){
               myIds <- which(detector.sp@data[ ,i] == myValues[j])
               plot(detector.sp[myIds,], col = myCols[j], pch=19, cex=2, add=TRUE)
            }
            
            plot(AC.sp, add=TRUE)
            plot(AC.sp[detected,], col = "red", pch = 19, add = TRUE)
            
            lapply(which(detected), function(x){
               this.row <- y.all[x, ]
               this.det <- detector.sp[this.row>0, ]
               if(length(this.det)>0){
                  segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                            , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                            , col = "red")
               }})
         }
      }
      
      ### ---- 3.Plot individual cov effects on sigma -----
      if(!is.null(id.covs.sigma) & type != "BinomFromPoisson" & type != "Recovery")
      {
         AC.sp@data <- id.covs.sigma
         for(i in 1:dim(id.covs.sigma)[2])
         {
            myCov <- id.covs.sigma[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(heat.colors(n = length(myValues), alpha = 1))
            
            par(mfrow=c(1,2))
            
            ## Plot the effect of covariate i on sigma
            myBeta <- id.betas.sigma[i]
            myX <- seq(min(myValues), max(myValues), length.out = 1000)
            if(link.sigma=="exponential"){myY <- exp(log(mySIGMA) + myBeta*myX)}
            else{myY <- mySIGMA + myBeta*myX}
            plot(myX, myY, xlab = names(id.covs.sigma[i]), ylab = "sigma", ylim = c(0,max(myY)), col = rev(heat.colors(length(myX))))
            
            ## Plot the Detections
            plot(AC.sp, main = names(myCov))
            plot(detector.sp, col="gray60", pch=19, cex=0.6, add=TRUE)
            
            for(j in 1:length(myValues)){
               myIds <- which(AC.sp@data[ ,i] == myValues[j])
               plot(AC.sp[myIds, ], col = myCols[j], pch = 19, add = TRUE)
               
               lapply(myIds,function(x){
                  this.row <- y.all[x, ]
                  this.det <- detector.sp[this.row>0, ]
                  if(length(this.det)>0){
                     segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                               , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                               , col = myCols[j])}})
            }
         }
      }
      
      ### ---- 4.Plot individual detections (no covariates) ----
      if(is.null(id.covs.sigma) & is.null(id.covs.p0) & is.null(det.covs.p0) & type != "Recovery" | type == "BinomFromPoisson")
         {
         if(!is.null(habitat.poly)){plot(habitat.poly)}else{plot(AC.sp)}
         plot(AC.sp,pch=1,cex=2, add=TRUE)
         plot(detector.sp, col="gray60", pch=19, cex=0.6, add = TRUE)
         plot(AC.sp[detected,], col = "red", pch = 19, add = TRUE)
         
         lapply(which(detected), function(x){
            this.row <- y.all[x, ]
            this.det <- detector.sp[this.row>0, ]
            if(length(this.det)>0){
               segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                         , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                         , col = "red")
            }})
         }
      
      ### ---- 5.Plot individual dead recoveries ----
      if(type == "Recovery")
         {
         plot(AC.sp)
         plot(detector.sp, col="gray60", pch=19, cex=0.6, add = TRUE)
         recovered <- which(y.all != max(y.all))
         
         plot(AC.sp[recovered,], col = "red", pch = 19, add = TRUE)
         lapply(recovered, function(x){
            this.det <- detector.sp[y.all[x], ]
            segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                    , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                    , col = "red")})
         }
      }
   ##------------------------------------------------------------------------------------------------------------------- 
   return(out)
   }
