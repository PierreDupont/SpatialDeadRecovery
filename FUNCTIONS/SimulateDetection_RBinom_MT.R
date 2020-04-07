
#' @title Function to simulate detection of simulated individuals within an SCR framework, under the influence of spatial covariates (at the detector level).
#'
#' @description
#' \code{SimulateDetection} returns a list object with \code{} 
#' 
#' @param sigma Numeric variable denoting the scale parameter of the detection function.
#' @param AC.sp Spatial points dataframe with individual activity center locations.
#' @param detector.sp Spatial points dataframe with detector locations.
#' @param detector.covs A data frame with the covariate values associated with each habitat grid cell (order identical to "covs").
#' @param betas A numeric vector denoting the coefficient values associated with the covariates to be used during simulations (effect on detection probability).
#' @param alpha A numeric variable denoting the probability that a sample is designated to an individual (i.e. implementation of partial identification)
#' @param n.samples A numeric variable denoting the total number of samples to be returned (only if \code{type}) is "Poisson". Output \code{p0} and \code{betas} are adjusted accordingly.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' areas on.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SimulateDetection.R
#' @keywords simul
#'
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' det.sim<-SimulateDetection(p0=0.6,sigma=1200,AC.sp=AC.sp,detector.sp=detector.sp,covs=covs, betas=betas,plot=TRUE)

SimulateDetection_MT<-function( p0 
                              , sigma
                              , AC.sp
                              , individual.covs
                              , individual.betas
                              , detector.sp
                              , detector.covs
                              , detector.betas 
                              , n.samples=NULL
                              , alpha=NULL
                              , type="RealBinom"
                              , plot=TRUE
                              , seed=NULL)
{
  X <- NULL
  projection(AC.sp)<-projection(detector.sp)  #---Just in case...
  
  #---occasionally rownames pose a problem; fix now:
  dimnames(detector.sp@data)[[1]]<-dimnames(detector.sp@coords)[[1]]<-1:length(detector.sp)
  dimnames(AC.sp@data)[[1]]<-dimnames(AC.sp@coords)[[1]]<-1:length(AC.sp)
  
  #-- Calculate the Detectors-Acs distance matrix
  D <- gDistance(detector.sp,AC.sp,byid=TRUE)
  
  #-- Calculate the logit of the intercept
  fixed.effects <- rep(logit(p0),length(detector.sp))
  
  #-- Calculate the effects for detector-specific covariates
  if(is.null(detector.covs)){detector.fixed.effects <- 0}
  if(!is.null(detector.covs)){
    temp <- formula(paste("~", paste(names(detector.covs), collapse="+"), sep=""))
    Xmat <- model.matrix(temp, detector.covs)
    detector.fixed.effects <- Xmat[,] %*% detector.betas
    detector.fixed.effects <- do.call(rbind, lapply(1:length(AC.sp), function(x) as.vector(t(detector.fixed.effects))))
    X.det <- data.frame(Xmat)
    }
  
  #-- Calculate the effects for individual covariates
  if(is.null(individual.covs)){individual.fixed.effects <- 0}
  if(!is.null(individual.covs)){
    temp <- formula(paste("~",paste(names(individual.covs),collapse="+"),"-1",sep=""))
    Xmat <- model.matrix(temp, individual.covs)
    individual.fixed.effects <- Xmat%*%individual.betas
    individual.fixed.effects <- do.call(cbind,lapply(1:length(detector.sp), function(x) as.vector(t(individual.fixed.effects))))
    X.ind <- data.frame(Xmat)
    }

  #-- Calculate the Baseline detection probability (including effects of detectors & individual covariates)
  P0 <- inv.logit(fixed.effects + individual.fixed.effects + detector.fixed.effects)
    
  #-- Calculate the sex-specific sigma
  if(length(sigma)>1){
      sigmaa <- D
      for(j in 1:nrow(D)){
        sigmaa[j, ] <- rep(sigma[(individual.covs[j,1]+1)], ncol(D)) # the detectors are 0 or one, so +1 is to avoid indexing  on zero.
         }
      sigma <-sigmaa
      }
  
  #-- Calculate the Detector-specific detection probability 
  P <- P0*exp(-D*D/(2*sigma*sigma))

  #-- Sample the individual & detector specific number of successes/detections
  y <- P
  for( i in 1:nrow(P))
     {
     for(j in 1:ncol(P))
        {
        y[i,j] <- rbinom(1,detector.sp$n.trials[j], P[i,j])
        }
     }
    
   dimnames(y)<-list(1:length(AC.sp),1:length(detector.sp))  #---DO NOT REMOVE! :)

  if(plot){
    AC.sp@data <- data.frame(sex=individual.covs )
    par(mfrow=c(1,2))
    
    AC1 <- AC.sp[AC.sp$sex==1,]
    y1 <- y[AC.sp$sex==1,]
    plot(AC1, col="pink", pch=19)
    plot(detector.sp,col="blue",pch=19,cex=0.6,add=TRUE)
    x<-1
    lapply(1:length(AC1),function(x){
       print(x)
       this.row <- y1[x,]
       this.det <- detector.sp[this.row>0, ]
       if(length(this.det)>0){
          segments(coordinates(this.det)[,1],coordinates(this.det)[,2],coordinates(AC1[x,])[,"x"],coordinates(AC1[x,])[,"y"],col="pink")
          }
    })
    
    AC0 <- AC.sp[AC.sp$sex==0,]
    y0 <- y[AC.sp$sex==0,]
    
    plot(AC0, col="lightblue", pch=19)
    plot(detector.sp,col="blue",pch=19,cex=0.6,add=TRUE)
    x<-1
    lapply(1:length(AC0),function(x){
       print(x)
       this.row<-y0[x,]
       this.det<-detector.sp[this.row>0,]
       if(length(this.det)>0){
          segments(coordinates(this.det)[,1],coordinates(this.det)[,2],coordinates(AC0[x,])[,"x"],coordinates(AC0[x,])[,"y"],col="lightblue")
         }
      })
   }
  
  #----SLIM TO INDIVIDUALS DETECTED
  y.all <- y
  detected <- apply(y,1,max)>0
  y <- y[detected,]
  
  out<-list(
    y=y,   #---this is the y for analysis
    y.all=y.all,#---this is the y to match with simulated ACs (contains all individuals, including those that were never detected)
    D=D,
    p0=p0, #lambda0 = lambda0
    sigma=sigma,
    betas.ind = individual.betas,
    betas.det = detector.betas,
    detector.covs=detector.covs,
    individual.covs=individual.covs,
    #det.sp=det.sp,
    alpha=alpha, # probability of making an ID
    n.trials = n.trials
  )  
  if(!is.null(detector.covs))out$X.det=X.det
  if(!is.null(individual.covs))out$X.ind=X.ind
  
  return(out)
}

