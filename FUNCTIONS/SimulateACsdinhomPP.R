#' @title Function to simulate activity center locations within a habitat mask, under the influence of spatial covariates using rinhomPP
#'
#' @description
#' \code{SimulateACs} returns a spatial point data frame with \code{N} simulated activity center 
#' locations within the space defined by \code{habitat.sp} and as function of the effects \code{betas} of spatial covariates rasters \code{cov}.
#' it also returns a matrix with the columns corresponding to the different rasters. A column (a vector) can be used as an input for the covariates in Nimble. 
#' 
#' @param N An \code{integer} denoting the number of simulated activity centers (ACs).
#' @param habitat.list object containing habitat list.
#' @param covs A \code{RasterLayer} object or rasterstack object when several covariates are to be used.
#' @param betas A numeric vector denoting the coefficient values associated with each raster.
#' @param seed if results to be reproduced, use a \code{numeric} value for the seed.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 
#'
#' @examples
#' # Generate a spatial point data frame with simulated activity center locations
#' mySimulatedACs <- SimulateACsdinhomPP(N= mySim$POPULATION$N
#'                                      , habitat.list= myHabitat.list
#'                                      , covs = MyDensityRaster
#'                                      , betas = mySim$DENSITY$Beta.dens
#'                                      , plot = TRUE
#'                                      , seed = NULL)

SimulateACsdinhomPP <- function(  N
                        , covs = NULL
                        , betas = NULL
                        , habitat.list = NULL
                        , plot = TRUE
                        , seed = NULL)
{
  
  habQualityDims <- dim(covs)[2:1]
  habQualityext <- as.vector(extent(covs))
  habQualityValues.arr <- array(NA, dim = dim(covs) )
  habQualityCovariate.arr <- array(NA, dim = dim(covs) )
  for( i in 1:dim(covs)[3]){
  # Convert raster to a matrix 
    habQualityCovariate.arr[,,i] <- as.matrix(covs[[i]])
  
  # Estimate habitat quality values from the regression coefficients
  habQualityValues.arr[,,i] <- matrix(
    exp(as.double(habQualityCovariate.arr[,,i]) * betas[i]),
    ncol = habQualityDims[1], nrow = habQualityDims[2])
  } 
  
  habQualityValues <- apply(habQualityValues.arr,c(1,2), sum)

  ## Set non-habitat cells as a 0 habitat quality value to avoid AC placement there.
  habQualityValues*habitat.list$habitat.mx
  
  # rescale the quality values so that the first covariates starts bottomright and fits dinhomPP
  habQualityValues <- habQualityValues[nrow(habQualityValues):1,]
 
  # Rearrange the habitat quality values to a vector so that the x-axis varies the fastest
  rHabQualityValues <- intensityRearrange(habQualityValues, habQualityDims, c(2, 1), c(0, 0))
  rHabQualityValues[is.na(rHabQualityValues)] <- 0

  # Simulate a set of activity centres drawn according to the habitat quality surface
  if(!is.null(seed)) {set.seed(seed+i)}
  AC <- t(replicate(N,
          # Use the inhomogenous Poisson-process to generate the activity centre (defined as a custom distribution)
          # IMPORTANT: Note here the use of the rearranged habitat quality values.
          rinhomPP(1
                   , intensityValues = rHabQualityValues
                   , intensityDims = habQualityDims
                   , lowerCoords = c(habQualityext[1], habQualityext[3])
                   , upperCoords = c(habQualityext[2], habQualityext[4]))
  ))
  
  # Convert to a sp object
  colnames(AC) <- c("x", "y")
  AC.sp <- data.frame(AC)
  coordinates(AC.sp) <- AC.sp
  AC.sp$Id <- 1:N
  
  ## plot the function 
  if(plot){
    # Plot the AC
    par(mfrow=c(1,2), mar=c(3,3,1,5))
    for( i in 1:dim(covs)[3]){
    plot(covs[[i]])
    points(AC.sp)
    
    # Plot the cov effect
    par( mar=c(5,5,1,1))
    
    myX <- seq(min(covs[[i]][],na.rm = T), max(covs[[i]][],na.rm = T), length.out = 1000)
    myY <- exp(betas[i]*myX)
    plot(myX, myY, xlab = names(covs)[i], ylab="intensity values", ylim = c(0,max(myY)), pch=19, col = rev(heat.colors(length(myX))))
    }

  }
  
  # Rearrange the habitat quality covariate to be used in the model
  HabQualityCovariate <- matrix(NA,prod(habQualityDims),dim(covs)[3])
  for( i in 1:dim(covs)[3]){
  HabQualityCovariate[,i] <- intensityRearrange(habQualityCovariate.arr[,,i], habQualityDims, c(2, 1), c(0, 0))
  }
  
  # create habitat toggle
  habitat.toggle <- as.numeric(!is.na(HabQualityCovariate[,1]))
  
  # replace NA by 0
  HabQualityCovariate[is.na(HabQualityCovariate)] <- 0
  
  return(list(   AC.sp = AC.sp
               , HabQualityCovariate = HabQualityCovariate
               , habitat.toggle = habitat.toggle
               )
         )
  
}

