### 10.1 ==== Define the density function ====
dbinomTeleportZ <- nimbleFunction(
  run = function(
    x = double(1),                               # Coordinate values to calculate the density
    lowerCoords = double(2),                     # The lower coordinate values of the observation windows
    upperCoords = double(2),                     # The upper coordinate values of the observation windows
    sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
    normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
    intensityWeightsMove = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    intensityWeightsRecruit = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    intensityWeightsTeleport = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
    numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
    dispersalToggle = double(0),    # A toggle (0/1) If individual performs large range dispersal or not
    firstAlive = double(0),
    Alive = double(0),
    localEvalParam = double(0, default = -1),    # A numeric specifying the distance to be use for the LESS approach
    log = integer(0, default = 0)                # If not 0 then return the log density
  ) {
    
    returnType(double(0))
    

    ## IF INACTIVE
    if(Alive == 0){
      # if(log){return(0)}else{return(1)}
       temp <- dbinomPPSingle(x,lowerCoords, upperCoords, intensityWeightsTeleport, areAreas, numWindows,log)
       return(temp)
    }

    ## FIRST TIME ALIVE
    if(firstAlive == 1){
      temp <- dbinomPPSingle(x,lowerCoords, upperCoords, intensityWeightsRecruit, areAreas, numWindows,log)
      return(temp)
    }

     if(Alive == 1){
    ## SECOND TIME ALIVE OR MORE
       if(dispersalToggle == 0){
         temp <- dbinomMNormSourcePPSingle(x, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeightsMove, areAreas, numWindows, localEvalParam, log)
         return(temp)
       }else{## TELEPORTATION 
         temp <- dbinomPPSingle(x,lowerCoords, upperCoords, intensityWeightsTeleport, areAreas, numWindows,log)
         return(temp)
       }
     }
    
  }
  
  
  
)



### 10.2. ==== Define the sampling function ====
rbinomTeleportZ <- nimbleFunction(
  run = function(
    n = integer(0),                              # Number of samples to draw from the distribution
    lowerCoords = double(2),                     # The lower coordinate values of the observation windows
    upperCoords = double(2),                     # The upper coordinate values of the observation windows
    sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
    normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
    intensityWeightsMove = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    intensityWeightsRecruit = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    intensityWeightsTeleport = double(1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
    areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
    numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
    dispersalToggle = double(0),
    firstAlive = double(0), 
    Alive = double(0),
    localEvalParam = double(0, default = -1)
  ) {
    returnType(double(1))
    if(dispersalToggle == 0){
      return(rbinomMNormSourcePPSingle(1, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeightsMove, areAreas, numWindows, localEvalParam)) 
    }else{
      return(rbinomPPSingle(1,lowerCoords, upperCoords, intensityWeightsMove, areAreas, numWindows)) 
      
    }
  }
)
##-------------------------------------------------------------------------------------------------------------------------------
## ------ XI.REGISTER NIMBLE DISTRIBUTIONS ------
### 11.1. ==== Register dbinomPP ====

registerDistributions(list(dbinomTeleportZ = list(
  BUGSdist = "dbinomTeleportZ(lowerCoords, upperCoords, sourceCoords, normSD, intensityWeightsMove, intensityWeightsRecruit, intensityWeightsTeleport , areAreas, numWindows, dispersalToggle, Alive, firstAlive, localEvalParam)",
  types = c(
    "value = double(1)", "lowerCoords = double(2)",
    "upperCoords = double(2)", "sourceCoords = double(1)", "normSD = double(0)",
    "intensityWeightsMove = double(1)", "intensityWeightsRecruit = double(1)","intensityWeightsTeleport = double(1)",
    "areAreas = double(0)", "numWindows = double(0)",
    "dispersalToggle = double(0)", "Alive = double(0)", "firstAlive = double(0)", "localEvalParam = double(0)" ),
  pqAvail = FALSE
)))
