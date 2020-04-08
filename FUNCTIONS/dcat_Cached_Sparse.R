#### 1.Density function ####
dcat_Cached_Sparse <- nimbleFunction(run = function( x = double(0),
                                                     pZero = double(0),
                                                     sigma = double(0),
                                                     sxy = double(1),
                                                     detectorCoords = double(2),
                                                     detectorID = double(2),
                                                     detectorNum = double(1),
                                                     habitatID = double(2), 
                                                     habitatFactor = double(0, default = 1.0),
                                                     indicator = double(0, default = 1.0),
                                                     log = double(0, default = 0.0)){
   ## RETURN TYPE DECLARATION
   returnType(double(0))
   
   ## GET NECESSARY INFO 
   n.detectors <- dim(detectorCoords)[1]
   alpha <- -1.0 / (2.0 * sigma * sigma)
   
   ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
   if(indicator == 0){
      if(x == 0){
         if(log == 0) return(1.0)
         else return(0.0)
      } else {
         if(log == 0) return(0.0)
         else return(-Inf)
      }
   }
   
   ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
   sID <- habitatID[trunc(sxy[2]/habitatFactor)+1, trunc(sxy[1]/habitatFactor)+1]
   
   ## GET NUMBER OF DETECTORS WITHIN maxDist FROM THE DETECTOR NUMBER MATRIX
   detNum <- detectorNum[sID]
   
   ## GET IDs OF DETECTORS WITHIN maxDist OF THE HABITAT CELL FROM THE DETECTOR ID MATRIX
   detIDs <- detectorID[sID, 1:detNum]
   
   ## Check that the recovery occurs at one of the allowed detectors
   ## if recovery outside allowed detectors logProb = -Inf (p = 0)
   if(sum(x == detIDs) == 0){
      if(log == 0) return(0.0)
      else return(-Inf)
   }
   
   ## CALCULATE THE LOG-LIKELIHOOD OF THE DETECTION
   ## Initialize the vector of detection probabilities (p = 0)
   p <- nimNumeric(length = detNum, value = 0.0, init = TRUE)
   newCat <- 0
   ## Loop over relevant detectors only
   for(j in 1:detNum){
      ## Calculate squared distance
      d2 <- pow(detectorCoords[detIDs[j],1] - sxy[1], 2) + pow(detectorCoords[detIDs[j],2] - sxy[2], 2)
      ## Calculate detection probability
      p[j] <- pZero * exp(alpha * d2)
      ## Establish the position of the detector in p[1:detNum]
      newCat <- newCat + (j * (x == detIDs[j]))
   }#j
   
   # sumP <- sum(p[1:detNum])
   # Calculate the log-likelihood of category x over n.categories 
   logProb <- dcat(newCat, prob = p[1:detNum], log = TRUE)
   if(log)return(logProb)
   return(exp(logProb))
   
   # logProb <- log(p[x]) - log(sumP)
   # if(log)return(logProb)
   # return(exp(logProb))
})

#### 2.Sampling function ####
rcat_Cached_Sparse <- nimbleFunction(run = function( n = double(0),
                                                     pZero = double(0),
                                                     sigma = double(0),
                                                     sxy = double(1),
                                                     detectorCoords = double(2),
                                                     detectorID = double(2),
                                                     detectorNum = double(1),
                                                     habitatID = double(2), 
                                                     habitatFactor = double(0, default = 1.0),
                                                     indicator = double(0, default = 1.0)){
   ## RETURN TYPE DECLARATION
   returnType(double(0))
   
   ## GET NECESSARY INFO 
   n.detectors <- dim(detectorCoords)[1]
   alpha <- -1.0 / (2.0 * sigma * sigma)
   
   ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
   if(indicator == 0){return(0)}
   
   ## GET HABITAT CELL ID FROM THE HABITAT ID MATRIX
   sID <- habitatID[trunc(sxy[2]/habitatFactor)+1, trunc(sxy[1]/habitatFactor)+1]
   
   ## GET NUMBER OF DETECTORS WITHIN maxDist FROM THE DETECTOR NUMBER MATRIX
   detNum <- detectorNum[sID]
   
   ## GET IDs OF DETECTORS WITHIN maxDist OF THE HABITAT CELL FROM THE DETECTOR ID MATRIX
   detIDs <- detectorID[sID, 1:detNum]
   
   ## INITIALIZE THE VECTOR OF DETECTIONS
   p <- nimNumeric(length = n.detectors, value = 0.0, init = TRUE)
   ## Loop over the relevant detectors only
   for(j in 1:detNum){
      ## Calculate squared distance
      d2 <- pow(detectorCoords[detIDs[j],1] - sxy[1], 2) + pow(detectorCoords[detIDs[j],2] - sxy[2], 2)
      ## Calculate detection probability
      p[detIDs[j]]  <- pZero * exp(alpha * d2)
   }#j
   
   ## OUTPUT
   detectOut <- rcat(n = 1, prob = p[1:n.detectors])
   return(detectOut)
})

#### 3.Registration ####
registerDistributions(list(
   dcat_Cached_Sparse = list(
      BUGSdist = "dcat_Cached_Sparse(pZero, sigma, sxy,
      detectorCoords, detectorID, detectorNum,
      habitatID, habitatFactor, indicator)",
      types = c( "value = double(0)","pZero = double(0)", "sigma = double(0)", "sxy = double(1)",
                 "detectorCoords = double(2)", "detectorID = double(2)", "detectorNum = double(1)",
                 "habitatID = double(2)", "habitatFactor = double(0)", "indicator = double(0)"),
      pqAvail = FALSE,
      mixedSizes = TRUE )))
