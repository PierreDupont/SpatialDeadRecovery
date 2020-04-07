#' @title Function to create a NIMBLE custom distribution for faster SCR model runs. i
#'
#' @description
#' \code{dbern_LESS} returns the likelihood of a given spatial detection history y[i,1:n.detectors] 
#'
#' @param x \code{Numeric} a scalar with the ID of the detector where detections happened (dcat).
#' @param pZero \code{Numeric} denoting the baseline detection probability of the half-normal at each detector.
#' @param sigma \code{Numeric} scalar denoting the scale parameter of the half-normal.
#' @param d2 \code{Vector} denoting the distance to each detector.
#' @param maxDist \code{Numeric} with the maximum distance to the AC where detections are considered possible.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{Integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dbern_LESS_vector(p0[1:n.detectors], sigma, d2[i,1:n.detectors,t], 10, z[i,t]==2)

#### 1.Density functions ####
dbern_LESS_singleCov <- nimbleFunction(run = function( x = double(1),
                                                       sxy = double(1),
                                                       pZero = double(0),
                                                       sigma = double(0),
                                                       detCoords = double(2),
                                                       detCov = double(1),
                                                       detBeta = double(0),
                                                       maxDist = double(0, default = 0.0),
                                                       indicator = double(0, default = 1.0),
                                                       log = integer(0, default = 0)){
   ## Return type declaration
   returnType(double(0))
   
   ## Check input dimensions
   n.detectors <- dim(detCoords)[1]
   
   ## Shortcut if individual is not available for detection
   if(indicator == 0){
      if(sum(x[1:n.detectors]) == 0){
         if(log == 0) return(1.0)
         else return(0.0)
      } else {
         if(log == 0) return(0.0)
         else return(-Inf)
      }
   }
   
   ## Process input values
   alpha <- -1.0 / (2.0 * sigma * sigma)
   maxDist_squared <- maxDist * maxDist
   
   ## Calculate detector-specific detection probabilities (only within maxDist of the AC ; else p = 0)
   p <- numeric(length = n.detectors, value = 0)
   for(j in 1:n.detectors){
      d2 <- pow(detCoords[j,1] - sxy[1], 2) + pow(detCoords[j,2] - sxy[2], 2)
      if(d2 <= maxDist_squared){ 
      p[j] <- exp(alpha * d2)/(1 + exp(-(logit(pZero) + detBeta * detCov[j])))
      }
   }#j
   
   ## Calculate the log-ikelihood of the given detection history
   logProb <- sum(dbinom(x, prob = p, size = 1, log = TRUE))
   if(log)return(logProb)
   return(exp(logProb))
})

#### 2.Sampling functions ####
rbern_LESS_singleCov <- nimbleFunction(run = function( n = integer(0),
                                                       sxy = double(1),
                                                       pZero = double(0),
                                                       sigma = double(0),
                                                       detCoords = double(2),
                                                       detCov = double(1),
                                                       detBeta = double(0),  
                                                       maxDist = double(0, default = 0.0),
                                                       indicator = double(0, default = 1.0)){
   ## Return type declaration
   returnType(double(1))
   
   ## Check input dimensions
   n.detectors <- dim(detCoords)[1]
   
   ## Shortcut if individual is not available for detection
   if(indicator == 0){return(rep(0.0, n.detectors))}
   
   ## Process input values
   alpha <- -1.0 / (2.0 * sigma * sigma)
   maxDist_squared <- maxDist * maxDist
   
   ## Calculate detector-specific detection probabilities (only within maxDist of the AC ; else p = 0)
   p <- numeric(length = n.detectors, value = 0)
   for(j in 1:n.detectors){
      d2 <- pow(detCoords[j,1] - sxy[1], 2) + pow(detCoords[j,2] - sxy[2], 2)
      if(d2 <= maxDist_squared){ 
         p[j] <- exp(alpha * d2)/(1 + exp(-(logit(pZero) + detBeta * detCov[j])))
      }
   }#j
   
   ## Draw detections from a Bernoulli distribution with the calculated probability
   y <- rbinom(n.detectors, 1, p)
   
   ## Output
   return(y)
})

#### 3.Registration ####
registerDistributions(list(
   dbern_LESS_singleCov = list(
      BUGSdist = "dbern_LESS_singleCov(sxy, pZero, sigma, detCoords, detCov, detBeta, maxDist, indicator)",
      types = c("value = double(1)",
                "sxy = double(1)",
                "pZero = double(0)",
                "sigma = double(0)",
                "detCoords = double(2)",
                "detCov = double(1)",
                "detBeta = double(0)",
                "maxDist = double(0)",
                "indicator = double(0)"),
      pqAvail = FALSE)))