### 10.1 ==== Define the density function ====
dbinomTeleport3 <- nimbleFunction( run = function( x = double(1),                               # Coordinate values to calculate the density
                                                   lowerCoords = double(2),                     # The lower coordinate values of the observation windows
                                                   upperCoords = double(2),                     # The upper coordinate values of the observation windows
                                                   sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
                                                   normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
                                                   intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
                                                   areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
                                                   numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
                                                   localEvalParam = double(0, default = -1),    # A numeric specifying the distance to be use for the LESS approach
                                                   indicator = double(0, default = 1),
                                                   log = integer(0, default = 0)){
  returnType(double(0))
  if(indicator == 0){
    # if(log){return(0)}else{return(1)}
    temp <- dbinomPPSingle(x,lowerCoords, upperCoords, intensityWeights, areAreas, numWindows,log)
    
  }else{
      temp <- dbinomMNormSourcePPSingle(x, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam, log)
      return(temp)
    
    }
})

### 10.2. ==== Define the sampling function ====
rbinomTeleport3 <- nimbleFunction( run = function(
  n = integer(0),                              # Number of samples to draw from the distribution
  lowerCoords = double(2),                     # The lower coordinate values of the observation windows
  upperCoords = double(2),                     # The upper coordinate values of the observation windows
  sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
  normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
  intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
  areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
  numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
  localEvalParam = double(0, default = -1),
  indicator = double(0, default = 1)
) {
  returnType(double(1))
  if(indicator == 1){
    return(rbinomMNormSourcePPSingle(1, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam))
  }else{
    return(rbinomPPSingle(1,lowerCoords, upperCoords, intensityWeights, areAreas, numWindows))
    
  }
}
)
##-------------------------------------------------------------------------------------------------------------------------------
## ------ XI.REGISTER NIMBLE DISTRIBUTIONS ------
### 11.1. ==== Register dbinomPP ====
registerDistributions(list(dbinomTeleport3 = list(
  BUGSdist = "dbinomTeleport3(lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam, indicator)",
  types = c(
    "value = double(1)", "lowerCoords = double(2)",
    "upperCoords = double(2)", "sourceCoords = double(1)", "normSD = double(0)", 
    "intensityWeights = double(1)", "areAreas = double(0)", "numWindows = double(0)", "localEvalParam = double(0)", "indicator = double(0)" ),
  pqAvail = FALSE
)))