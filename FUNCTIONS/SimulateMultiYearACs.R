#' @title Function to simulate several individual AC locations within a habitat mask over multiple years
#'
#' @description
#' \code{SimulateMultiYearACs} returns a list of length\code{n.occasions} containing \code{SpatialPointDataFrame} with \code{N} simulated AC locations within the space defined by \code{habitat.sp}. 
#' requires package {circular}
#' 
#' @param N A \code{integer}: number of individuals for which to simulate ACs.
#' @param n.occasions A \code{integer}: number of years (=time steps) for which to simulate ACs.
#' @param habitat.sp A \code{SpatialPointDataFrame} specifying the habitat mask (grid centers).
#' @param type.SS A \code{character} chain specifying the nature of the State-Space to be used. Can be one of \code{"Continuous"} or \code{"Discrete"}
#' @param type A \code{character} chain specifying the method to be used. Can be one of \code{"Hyper.AC"}; \code{"Markov.Angle"}; \code{"Markov.XY"} or \code{"Markov.cov}
#'
#' HYPER.AC OPTION
#' @param sigma A \code{numeric} only used if type == "Hyper.AC" ; standard deviation of the distance between each year AC and the hyper-AC location.
#'
#' MARKOV.COV OPTION
#' @param sigma.hat A \code{numeric} only used if type == "Markov.cov" ; standard deviation of the half-normal movement distribution.
#' 
#' MARKOV.ANGLE OPTION
#' @param mean.angle A \code{numeric} only used if type == "Markov.Angle"; mean angle (in radians) of the direction of the movement between two consecutive time steps.
#' @param sd.angle A \code{numeric} only used if type == "Markov.Angle"; standard deviation (in radians) of the direction of the movement between two consecutive time steps.
#' @param distance.distribution A \code{character} only used if type == "Markov.Angle"; sets the distribution to be used for sampling distances between two consecutive ACs; can be one of \code{"Poisson"}, \code{"Exponential"}, \code{"Log.normal"} or \code{"Gamma"}.
#' @param lambda.pois A \code{numeric} only used if distance.distribution == "Poisson"; mean of the Poisson distribution 
#' @param rate.exp A \code{numeric} only used if distance.distribution == "Exponential"; rate of the Exponential distribution 
#' @param rate.gamma A \code{numeric} only used if distance.distribution == "Gamma"; rate of the Gamma distribution 
#' @param shape.gamma A \code{numeric} only used if distance.distribution == "Gamma"; shape parameter of the Gamma distribution 
#' @param mean.lnorm A \code{numeric} only used if distance.distribution == "Log.normal"; mean of the Log.normal distribution 
#' @param sd.lnorm A \code{numeric} only used if distance.distribution == "Log.normal"; standard deviation of the Log.normal distribution 
#' 
#' MARKOV.XY OPTION
#' @param mean.norm.XY A \code{numeric} vector containing the mean movement distances (normal distributions) along X & Y respectively. If only one value is provided it applies to both X and Y distance sampling. 
#' @param sd.norm.XY A \code{numeric} vector containing the standard deviations of the movement distances (normal distributions) along X & Y. If only one value is provided it applies to both X and Y distance sampling. 
#' 
#' OPTIONAL PARAMETERS
#' @param plot.check A \code{logical} for whether \code{TRUE} or \code{FALSE} plots should be displayed to check simulations. 
#' @param initial.AC.sp A \code{SpatialPointDataFrame} specifying the initial individual AC locations.
#' @param covs A \code{DataFrame}: covariate values associated with each habitat grid cell (order identical to habitat.sp).
#' @param betas A \code{numeric} vector: regression coefficients associated with each \code{covs}
#' @param habitat.poly A \code{SpatialPolygon} of the habitat mask for plotting (and eventually initial locations sampling if type.SS is continuous)
#' 
#' @author Pierre Dupont, \email{pierre.dupont@@nmbu.no}
#' @backref R/SimulateMultiYearACs.R
#' @keywords simulation
#'
#' @examples
#' # Generate a habitat spatial point data frame with simulated activity center locations over multiple years.
#' tab <- data.frame("x"= sample(100,500,replace=TRUE), "y"=sample(200,500,replace=TRUE))
#' habitat.sp <- SpatialPointsDataFrame(tab, data = tab)
#'
#' # Acs normally distributed around a Hyper-AC:
#' AC.ls <- SimulateMultiYearACs( N = 20
#'                              , n.occasions = 10
#'                              , habitat.sp =  habitat.sp
#'                              , type = "Hyper.AC"
#'                              , sigma = 2
#'                              , plot.check = TRUE)
#'                              
#' AC.ls <- SimulateMultiYearACs( N = 20
#'                               , n.occasions = 10
#'                               , habitat.sp =  habitat.sp
#'                               , type = "Markov.Angle"
#'                               , mean.angle = 0
#'                               , sd.angle = 1
#'                               , distance.distribution = "Log.normal"
#'                               , mean.lnorm = 0.2
#'                               , sd.lnorm = 0.1)

SimulateMultiYearACs <- function( N 
                                , n.occasions
                                , habitat.sp
                                , type.SS = "Discrete"
                                , type = "Hyper.AC"
                                ## If type = "Hyper.Ac"
                                , sigma = NULL
                                ## If type = "Markov.Cov"
                                , sigma.hat = NULL
                                ## If type = "Markov.Angle"
                                , mean.angle = pi/2
                                , sd.angle = 0.1
                                , distance.distribution = "Poisson"
                                , lambda.pois = 1
                                , rate.exp = 1
                                , rate.gamma = 1
                                , shape.gamma = 1
                                , mean.lnorm = 1
                                , sd.lnorm = 1
                                ## If type = "Markov.XY"
                                , mean.norm.XY = c(1,1)
                                , sd.norm.XY = c(1,1)
                                ## OTHERs
                                , initial.AC.sp = NULL
                                , covs = NULL
                                , betas = NULL
                                , plot.check = TRUE
                                , habitat.poly = NULL)
   {
   ##-------------------------------------------------------------------------------------
   ## ==== STEP 1: Initialize AC positions ====
   if(is.null(habitat.poly)){
      min.x <- min(coordinates(habitat.sp)[,1])
      max.x <- max(coordinates(habitat.sp)[,1])
      min.y <- min(coordinates(habitat.sp)[,2])
      max.y <- max(coordinates(habitat.sp)[,2])
      
      habitat.poly <- MakePolygon(coord.x = c(min.x, min.x, max.x, max.x), coord.y = c(min.y, max.y, max.y, min.y))
      }#if
   
   AC.sp.list <- list()
      if(!is.null(initial.AC.sp)){ 
         N <- length(initial.AC.sp)
         AC.sp.list[[1]] <- initial.AC.sp
      }
 
   if(is.null(initial.AC.sp)){
      if(type.SS == "Discrete"){
         AC.sp.list[[1]] <- SimulateACs(N = N, habitat.sp = habitat.sp, covs = NULL, betas = NULL, plot = FALSE)
         AC.sp.list[[1]] <- AC.sp.list[[1]][sample(1:N,N),]
          }#if
      
      if(type.SS == "Continuous"){
         temp.AC <- spsample(x = habitat.poly, n = N, type = "random")
         AC.sp.list[[1]] <- SpatialPointsDataFrame( temp.AC
                                                  , data = data.frame(ID = 1:N)
                                                  , proj4string = CRS(projection(habitat.sp)))
         }#if

      rownames(AC.sp.list[[1]]@coords) <- 1:length(AC.sp.list[[1]])
      }#if
   

   ##-------------------------------------------------------------------------------------
   ################################################
   ### ------ OPTION 1 : Hyper-parameter ------ ###
   ################################################
   if(type == "Hyper.AC"){
      
      hyperAC.sp <- AC.sp.list[[1]]
      
      ## ==== 1.DISCRETE STATE-SPACE ====
      if(type.SS == "Discrete"){
         ## STEP 2: Calculate habitat-specific probabilities based on hyper-ACs 
         D <- gDistance(habitat.sp, hyperAC.sp, byid = TRUE)
         prob.D <- dnorm(x = D, mean = 0, sd = sigma)
         
         ## STEP 3: Generate ACs list
         AC.index <- apply(prob.D, 1, function(x)which(rmultinom(n.occasions, 1, x)==1, arr.ind = TRUE)[ ,1])
         AC.sp.list <- apply(matrix(AC.index, ncol = N, nrow = n.occasions ), 1, function(x){
                                                      coords <- coordinates(habitat.sp[x,])
                                                      row.names(coords) <- 1:N
                                                      SpatialPointsDataFrame( coords = coords
                                                                            , data = data.frame(ID = 1:N)
                                                                            , proj4string = CRS(projection(habitat.sp)))
                                           })
         
         }#if
      
      ## ==== 2.CONTINUOUS STATE-SPACE ==== 
      if(type.SS == "Continuous"){
         ## STEP 2: Sample individual dispersal distance along X & Y axes
         D.X <- sapply(1:n.occasions, function(x)rnorm(n = N, mean = 0, sd = sigma))
         D.Y <- sapply(1:n.occasions, function(x)rnorm(n = N, mean = 0, sd = sigma))

         ## STEP 3: Generate new AC position 
         myX <- coordinates(hyperAC.sp)[ ,1] + D.X
         myY <- coordinates(hyperAC.sp)[ ,2] + D.Y
         
         ## STEP 3: Store the new AC coordinates in a SpatialPointDataFrame 
         AC.sp.list <- lapply(1:dim(myX)[2], function(x)SpatialPointsDataFrame( coords = cbind.data.frame(x = myX[,x] , y = myY[,x]) 
                                                                              , data = data.frame(ID = 1:N)
                                                                              , proj4string = CRS(projection(habitat.sp))))
         }#if
      }#if
   
    ##-------------------------------------------------------------------------------------
   ######################################################################
   ### ------ OPTION 2 : Markovian Movement : ANGLE & DISTANCE ------ ###
   ######################################################################
   if(type == "Markov.Angle"){
      for(t in 2:n.occasions){
         ## ==== 1.DISCRETE STATE-SPACE ====
         if(type.SS == "Discrete"){
            ## STEP 2: Calculate distance and angles between each habitat cell and each AC
            D <- apply(coordinates(habitat.sp), 1, function(x){sqrt((coordinates(AC.sp.list[[t-1]])[ ,1]-x[1])^2 + (coordinates(AC.sp.list[[t-1]])[ ,2]-x[2])^2)})
            A <- apply(coordinates(habitat.sp), 1,function(x){atan2(x[2]-coordinates(AC.sp.list[[t-1]])[ ,2],x[1]- coordinates(AC.sp.list[[t-1]])[ ,1])})
            
            ## STEP 3: Generate individual & habitat-specific dispersal probabilities based on DISTANCE 
            if(distance.distribution == "Exponential"){prob.D <- dexp(x = D, rate = rate.exp)}
            if(distance.distribution == "Gamma"){prob.D <- dgamma(x = D, shape = shape.gamma, rate = rate.gamma)} 
            if(distance.distribution == "Log.normal"){prob.D <- dlnorm(x = D, meanlog = mean.lnorm, sdlog = sd.lnorm)}
            if(distance.distribution == "Poisson"){prob.D <- dpois(x = D, lambda = lambda.pois)}
            
            ## STEP 4: Generate individual & habitat-specific dispersal probabilities based on ANGLES
            suppressWarnings(prob.A <- dvonmises(x = A, mu = mean.angle, kappa = 1/(sd.angle^2), log = FALSE))
            
            ## STEP 5: Sample new AC position  
            AC.prob <- prob.A*prob.D
            if(!is.null(covs)){
               temp <- formula(paste("~", paste(names(covs), collapse="+"), sep=""))
               Xmat <- model.matrix(temp, covs)
               explambda <- exp(Xmat[,] %*% betas)
               mu <- explambda/sum(explambda) 
               }#if
            AC.index <- apply(AC.prob, 1, function(x)which(rmultinom(n = 1, size = 1, prob = x)==1))
            
            ## STEP 6: Store the new AC coordinates in a SpatialPointDataFrame 
            AC.sp.list[[t]] <- SpatialPointsDataFrame( coords = coordinates(habitat.sp[AC.index, ])
                                                     , data = data.frame(ID = 1:N)
                                                     , proj4string = CRS(projection(habitat.sp)))
            }#if
         
         ## ==== 2.CONTINUOUS STATE-SPACE ==== 
         if(type.SS == "Continuous"){
            ## STEP 2: Sample movement distances and angles ==== 
            if(distance.distribution == "Exponential"){D <- rexp(n = N, rate = rate.exp)}
            if(distance.distribution == "Gamma"){D <- rgamma(n = N, shape = shape.gamma, rate = rate.gamma)} 
            if(distance.distribution == "Log.normal"){D <- rlnorm(n = N, meanlog = mean.lnorm, sdlog = sd.lnorm)}
            if(distance.distribution == "Poisson"){D <- rpois(n = N, lambda = lambda.pois)}
            
            suppressWarnings(theta <- rvonmises(n = N, mu = mean.angle, kappa = 1/(sd.angle^2), control.circular = list(type = "angles", units = "radians", template = "none", modulo = "asis", zero = 0, rotation = "counter")))
   
            ## STEP 3: Generate new AC positions ==== 
            x <- coordinates(AC.sp.list[[t-1]])[ ,1] + sin(theta)*D
            y <- coordinates(AC.sp.list[[t-1]])[ ,2] + cos(theta)*D
            
            ## STEP 4: Store the new AC coordinates in a SpatialPointDataFrame 
            AC.sp.list[[t]] <- SpatialPointsDataFrame( coords = cbind.data.frame(x = x, y = y)
                                                     , data = data.frame(ID = 1:N)
                                                     , proj4string = CRS(projection(habitat.sp)))
            }#if
            
         rownames(AC.sp.list[[t]]@coords) <- 1:length(AC.sp.list[[t]])
         }#t
      }#if
   
   ##-------------------------------------------------------------------------------------
   #############################################################################
   ### ------ OPTION 3 : Markovian Movement : X.DISTANCE & Y.DISTANCE ------ ###
   #############################################################################
   if(type == "Markov.XY"){
      if(length(mean.norm.XY)==1){mean.norm.XY <- rep(mean.norm.XY,2)}
      if(length(sd.norm.XY)==1){sd.norm.XY <- rep(sd.norm.XY,2)}
      for(t in 2:n.occasions){
         ## ==== 1.DISCRETE STATE-SPACE ====
         if(type.SS == "Discrete"){
            ## STEP 2: Calculate distances along X & Y between each habitat cell and each AC 
            D.X <-  apply(coordinates(habitat.sp), 1, function(x){coordinates(AC.sp.list[[t-1]])[ ,1] - x[1]})
            D.Y <-  apply(coordinates(habitat.sp), 1, function(x){coordinates(AC.sp.list[[t-1]])[ ,2] - x[2]})

            ## STEP 3: Generate individual & habitat-specific dispersal probabilities based on X & Y distances
            prob.X <- dnorm(x = D.X, mean = mean.norm.XY[1], sd = sd.norm.XY[1])
            prob.Y <- dnorm(x = D.Y, mean = mean.norm.XY[2], sd = sd.norm.XY[2])
            
            ## STEP 4: Sample new AC position
            AC.prob <- prob.X*prob.Y
            AC.index <- sapply(1:dim(AC.prob)[1], function(x)which(rmultinom(n = 1, size = 1, prob = AC.prob[x, ])==1))
            
            ## STEP 5: Store the new AC coordinates in a SpatialPointDataFrame
            AC.sp.list[[t]] <- SpatialPointsDataFrame( coords = data.frame("x"= coordinates(habitat.sp[AC.index, ])[,1],"y"= coordinates(habitat.sp[AC.index, ])[,2])
                                                     , data = data.frame(ID = 1:N)
                                                     , proj4string = CRS(projection(habitat.sp)))
            }#if
         
         ## ==== 2.CONTINUOUS STATE-SPACE ====
         if(type.SS == "Continuous"){
            ## STEP 2: Sample individual dispersal distance along X & Y axes
            D.X <- rnorm(n = N, mean = mean.norm.XY[1], sd = sd.norm.XY[1])
            D.Y <- rnorm(n = N, mean = mean.norm.XY[2], sd = sd.norm.XY[2])
            
            ## STEP 3: Generate new AC position 
            x <- coordinates(AC.sp.list[[t-1]])[ ,1] + D.X
            y <- coordinates(AC.sp.list[[t-1]])[ ,2] + D.Y
            
            ## STEP 4: Store new AC coordinates in a SpatialPointDataFrame 
            AC.sp.list[[t]] <- SpatialPointsDataFrame( coords = cbind.data.frame(x = x, y = y) 
                                                     , data = data.frame(ID = 1:N)
                                                     , proj4string = CRS(projection(habitat.sp)))
            }#if
         
         rownames(AC.sp.list[[t]]@coords) <- 1:length(AC.sp.list[[t]])
         }#t
      }#if
   
   ##-------------------------------------------------------------------------------------
   ###############################################################################
   ### ------ OPTION 4 : Markovian Movement : Half-normal with covariates ---- ###
   ###############################################################################
   if(type == "Markov.Cov"){
      for(t in 2:n.occasions){
         ## STEP 2: Calculate distances along between each habitat cell and each AC 
         D <- apply(coordinates(habitat.sp), 1, function(x){sqrt((coordinates(AC.sp.list[[t-1]])[ ,1] - x[1])^2 + (coordinates(AC.sp.list[[t-1]])[ ,2] - x[2])^2)})
         
         ## STEP 3: Generate individual & habitat-specific dispersal probabilities based on distance only
         prob.D <- dnorm(x = D, mean = 0, sd = sigma.hat)
         
         ## STEP 4: Generate habitat-specific effects
         fixed.effects.p <- 0
         if(!is.null(covs)){
            temp <- formula(paste("~", paste(names(covs), collapse = "+"),"-1"), sep = "")
            Xmat <- model.matrix(temp, covs)
            fixed.effects.p <- exp(Xmat %*% betas)
            }#if
         fixed.effects.p <- fixed.effects.p/sum(fixed.effects.p)
         prob.C <- do.call(rbind,lapply(1:dim(D)[1], function(x)t(fixed.effects.p)))

         ## STEP 4: Sample new AC position
         AC.prob <- prob.D*prob.C
         AC.index <- sapply(1:dim(AC.prob)[1], function(x)which(rmultinom(n = 1, size = 1, prob = AC.prob[x, ])==1))
            
         ## STEP 5: Store the new AC coordinates in a SpatialPointDataFrame
         AC.sp.list[[t]] <- SpatialPointsDataFrame( coords = data.frame("x"= coordinates(habitat.sp[AC.index, ])[,1],"y"= coordinates(habitat.sp[AC.index, ])[,2])
                                                  , data = data.frame(ID = 1:N)
                                                  , proj4string = CRS(projection(habitat.sp)))
            
         rownames(AC.sp.list[[t]]@coords) <- 1:length(AC.sp.list[[t]])
         }#t
      }#if
   
   ##-------------------------------------------------------------------------------------
   ## ==== STEP 5: GENERATE NEW HABITAT POLYGON AROUND ACs WHEN NEEDED ==== 
   min.x <- min(c(unlist(lapply(AC.sp.list,function(x){coordinates(x)[,1]})), extent(habitat.poly)[1]))
   max.x <- max(c(unlist(lapply(AC.sp.list,function(x){coordinates(x)[,1]})), extent(habitat.poly)[2]))
   min.y <- min(c(unlist(lapply(AC.sp.list,function(x){coordinates(x)[,2]})), extent(habitat.poly)[3]))
   max.y <- max(c(unlist(lapply(AC.sp.list,function(x){coordinates(x)[,2]})), extent(habitat.poly)[4]))
      
   new.poly <- MakePolygon(coord.x=c(min.x, min.x, max.x, max.x), coord.y=c(min.y,max.y,max.y,min.y))
   
   ##-------------------------------------------------------------------------------------
   ## ==== STEP 6: Plot ACs ====
   if(plot.check){
      
      if(type == "Markov.Cov"){
         cov.r <- rasterFromXYZ(habitat.sp)
         cov.r[!is.na(cov.r)] <- fixed.effects.p
         plot(cov.r)
      }else{if(any(is.na(unlist(lapply(AC.sp.list, function(x)over(x = x, y = habitat.poly)))))){
         plot(new.poly, col="red")
         plot(habitat.sp, col = grey(0.7), pch=19, cex=1, add=TRUE)
      }else{plot(habitat.sp, col = grey(0.7), pch=19, cex=1)}#else
      }#else
      
      plot(habitat.poly, add = TRUE)
      col <- colors()[sample(max(100,N),N)]
      
      if(type == "Hyper.AC"){
         lapply(1:length(AC.sp.list), function(x){
            plot(AC.sp.list[[x]], add = TRUE, col = col, pch = 19, cex = 1)
            segments(coordinates(AC.sp.list[[x]])[,1], coordinates(AC.sp.list[[x]])[,2]
                     , coordinates(hyperAC.sp)[,1], coordinates(hyperAC.sp)[,2]
                     , col = col)
            })
         plot(hyperAC.sp, add=TRUE, bg=col, col="black", lwd=3, pch=21, cex=2)
      }else{
         lapply(2:length(AC.sp.list), function(x){         
            suppressWarnings(arrows( x0 = AC.sp.list[[x-1]]$x, y0 = AC.sp.list[[x-1]]$y
                  , x1 = AC.sp.list[[x]]$x, y1 = AC.sp.list[[x]]$y 
                  , col = col, lwd = 2, length = 0.07, angle = 30))
            })
         lapply(AC.sp.list, function(x){plot(x, add = TRUE, bg=col, col="black", pch=21, cex=1)})
         }#else
      
   }#if
   
   ##-------------------------------------------------------------------------------------
   ## ==== STEP 7: Output ====
   out <- list(AC.sp.list = AC.sp.list)
   if(type == "Hyper.AC"){out$hyperAC.sp = hyperAC.sp}
   if(any(is.na(unlist(lapply(AC.sp.list, function(x)over(x = x, y = habitat.poly)))))){out$new.habitat.poly = new.poly}
   return(out)
   
   ##-------------------------------------------------------------------------------------
   }



