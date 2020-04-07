#' @title Function to simulate activity center locations within a habitat mask, under the influence of spatial covariates
#'
#' @description
#' \code{SimulateACs} returns a spatial point data frame with \code{N} simulated activity center locations within the space defined by \code{habitat.sp} and as function of the effects \code{betas} of spatial covariates \code{cov}.
#'  The projection of \code{habitat.sp} will be retained for the output spatial points data frame.
#' 
#' @param N An integer denoting the number of simulated activity centers (ACs).
#' @param habitat.sp A spatial points data frame specifying the habitat mask (grid centers).
#' @param covs A data frame with the covariate values associated with each habitat grid cell (order identical to habitat.sp).
#' @param betas A numeric vector denoting the coefficient values associated with the covariates to be used during simulations (effect on the multinomial probability of AC placement across the available habitat cells).
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SimulateACs.R
#' @keywords simul
#'
#' @examples
#' # Generate a spatial point data frame with simulated activity center locations
#' AC.sp<-SimulateACs(N=5,habitat.sp=habitat.ls$habitat.sp,covs=covs,betas=betas,plot=TRUE)

SimulateACs<-function( N
                     , habitat.sp
                     , covs = NULL
                     , betas = NULL
                     , plot = TRUE
                     , seed = NULL)
   {
   ## ---- MAKE MODEL MATRIX AND IMPLEMENT 
   mu <- rep(1,length(habitat.sp))
   
   if(!is.null(covs)){
      temp <- formula(paste("~", paste(names(covs), collapse="+"), sep=""))
      Xmat <- model.matrix(temp, covs)
      explambda <- exp(Xmat[,] %*% betas)
      mu <- explambda/sum(explambda) 
      }#if
   
   if(!is.null(seed)) {set.seed(seed)}
   
   # temp <- as.numeric(rmultinom(1, N, mu))
   # index <- rep(1:length(temp), temp)
   index <- replicate(expr = which(rmultinom(1,1,mu)==1), n = N)
   
   AC.sp <- habitat.sp[index, ] ###--if we don't want cell center points, could instead sample points randomly within each given habitat raster cell

   if(plot){
      plot(habitat.sp, pch=19, col="grey")
      points(jitter(y, 2) ~ jitter(x, 2), coordinates(AC.sp), col="red", pch=19)
      }#if
   
   return(AC.sp)
   }

