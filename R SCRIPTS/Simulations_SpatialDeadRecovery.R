#############################################################################
####### ---- Open-Population Spatial Capture-Recapture-Recovery ---- ########
####### ------ DATA SIMULATION & NIMBLE MODEL FITTING SCRIPT ------- ########
#############################################################################
rm(list=ls())

## ------ IMPORT REQUIRED LIBRARIES ------
library(rgeos)                    # Import the geospatial analysis libraries
library(raster)                    # Import the raster spatial package
library(coda)                      # Import the MCMC diagnostic tools
library(nimble)                    # Import the NIMBLE subroutines
library(R.utils)                   # Import the utility functions (some of the API of this package is experimental)
library(abind)                     # Import the library for manipulating multidimensional arrays

## ------ SOURCE THE REQUIRED FUNCTIONS ------
source("FUNCTIONS/pointProcess.R")       ## NIMBLE functions for the point process model
source("FUNCTIONS/MakeZ.R")              ## R function to generate random individual states
source("FUNCTIONS/MakeSxyInits.R")       ## R function to generate random individual activity center locations
source("FUNCTIONS/dbern_LESS.R")         ## NIMBLE custom Bernoulli distribution for a vector of binary observations at multiple detectors
source("FUNCTIONS/dcat_LESS.R")          ## NIMBLE custom categorical distribution for dead recovery observations
source("FUNCTIONS/calculateDistance.R")  ## NIMBLE custom function to calculate the squared distance between an activity center and all detectors

## ---------------------------------------------------------------------------------
## ------ I.SET SIMULATION PARAMETERS -----
sim <- list( "tau" = 5,     ## Standard deviation of the bivariate normal movement distribution
             "N0" = 50,     ## Initial population size
             "n.years" = 5, ## Number of simulated years
             "rep" = 0.4,   ## Reproduction probability
             "fec" = 1,     ## Mean per breeder recruits
             "phi" = 0.6,   ## Survival probability
             "r" = 0.5,     ## Recovery probability
             "p0" = 0.15,   ## Baseline detection probability
             "sigma" = 2)   ## Scale parameter of the detection function

## ---------------------------------------------------------------------------------
## ------ II.DATA SIMULATION ------
### ====    1. GENERATE HABITAT CHARACTERISTICS ====
## Generate a spatial polygon of the study area
studyAreaPolygon <- SpatialPolygons(Srl = list( Polygons(srl = list(Polygon(cbind( c(0,0,30,30), c(0,30,30,0)))),
                                                         ID = "studyArea")),
                                    proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Buffer the study area
bufferedStudyArea <- gBuffer(spgeom = studyAreaPolygon, width = 2.5)

## Generate a raster of available habitat
habitat.r <- raster( x = extent(bufferedStudyArea),
                     resolution = 17.5,
                     crs = proj4string(bufferedStudyArea))
habitat.r <- rasterize(bufferedStudyArea, habitat.r)

## Retrieve lower and upper cell coordinates
numHabCells <- length(habitat.r)
lowerCoords <- coordinates(habitat.r) - 0.5*res(habitat.r)
upperCoords <- coordinates(habitat.r) + 0.5*res(habitat.r)

## Plots
plot(habitat.r)
plot(bufferedStudyArea, add = T)
plot(studyAreaPolygon, add = T, col = "gray60")

### ====    2. GENERATE DETECTORS ====
## Generate a raster of detector cells
detector.r <- raster( x = extent(studyAreaPolygon),
                      resolution = 1.5,
                      crs = proj4string(studyAreaPolygon))
detector.r <- rasterize(studyAreaPolygon, detector.r)

## Generate data frame of detectors coordinates 
detector.xy <- as.data.frame(coordinates(detector.r))
points(detector.xy[ ,1], detector.xy[ ,2], pch = 3, cex = 0.5)

### ====    3. GENERATE POPULATION DYNAMICS ====
phi <- sim$phi
r <- sim$r
rep <- sim$rep
fec <- sim$fec
N0 <- sim$N0
n.years <- sim$n.years

## Build state transition matrix (3 states; 1:alive, 2:recovered, 3:dead)
ST <-  matrix(c( phi, (1-phi)*r, (1-phi)*(1-r),
                 0  , 0, 1, 
                 0  , 0, 1), nrow = 3, byrow = TRUE)

## Initialize list of population composition 
POP <- list()
POP[[1]] <- data.frame(id = 1:N0, z = rep(1,N0))

## Initialize vector of population size 
N <- NULL
N[1] <- N0

## Loop demographic process over n.years
for(t in 2:n.years){
  ## Sample number of breeders
  R <- rbinom(n = 1, size = N[t-1], rep)
  
  ## Sample number of recruits
  B <- rpois(n = R, fec)
  N.new <- sum(B[])
  
  ## Sample individual survival/state transition
  Z <- NULL
  for(i in 1:(dim(POP[[t-1]])[1])){
    Z[i] <- which(rmultinom(1, 1, ST[POP[[t-1]]$z[i], ]) == 1)
  }#i

  ## New population composition
  Z <- c(Z, rep(1,N.new))
  POP[[t]] <- data.frame(id = 1:length(Z), z = Z)
  
  ## Population size
  N[t] <- sum(Z == 1)
}#t

## Generate z matrix based on population list
z <- matrix(NA, length(Z), n.years)
for(t in 1:n.years){
  z[1:dim(POP[[t]])[1],t] <- POP[[t]]$z
}#t

## Add 1 for state "available" 
## (4 states; 1: available, 2: alive, 3: recovered, 4: dead)
true.z <- z
true.z <- true.z + 1 
true.z[is.na(true.z)] <- 1

### ====    4. GENERATE INDIVIDUAL AC LOCATIONS ====
## Initiliaze the array of activity center coordinates
sCoords <- array(NA, c(dim(true.z)[1], 2, n.years))

## Sample original AC locations
sCoords[ , ,1] <- rbinomPP( n = 1,
                            numPoints = dim(true.z)[1],
                            lowerCoords = lowerCoords,
                            upperCoords = upperCoords,
                            intensityWeights = rep(1, numHabCells),
                            areAreas = 1,
                            numWindows = nrow(lowerCoords))

## Sample AC movements
for(t in 2:dim(true.z)[2]){
  sCoords[ , ,t] <- rbinomMNormSourcePPMulti( n = 1,
                                              lowerCoords = lowerCoords,
                                              upperCoords = upperCoords,
                                              sourceCoords = sCoords[ , ,t-1],
                                              normSD = sim$tau,
                                              intensityWeights = rep(1, numHabCells),
                                              areAreas = 1,
                                              numWindows = nrow(lowerCoords))
  
}#t
points(sCoords[,1,], sCoords[,2,], col = "red", pch = 19, cex = 0.5)

### ====    5. GENERATE ALIVE DETECTIONS : y.alive[i,j,t] ====
## List individuals alive (= available for detection)
z.alive <- apply(true.z, 2, function(x){which(x == 2)})

## Initialize array of detections (individuals * detectors * years)
y.alive <- array(0, c(dim(true.z)[1], dim(detector.xy)[1], n.years))

## Sample individual detections
for(t in 1:n.years){
  #--- Check that at least one individual available for detection
  if(length(z.alive[[t]]) > 0){
    for(i in z.alive[[t]]){
      #--- Calculate squared distance to all detectors
      D2 <- (detector.xy[ ,1] - sCoords[i,1,t])^2 + (detector.xy[ ,2] - sCoords[i,2,t])^2
      #--- Calculate detector-specific detection probability (half-normal detection function)
      P <- sim$p0 * exp(-D2/(2 * sim$sigma * sim$sigma))
      #-- Sample binary individual detections 
      y.alive[i, ,t] <- rbinom(n = length(P), size = 1, P)
      points(detector.xy[y.alive[i, ,t] == 1,1],
             detector.xy[y.alive[i, ,t] == 1,2],
             pch = 3, col = "blue", cex = 0.5)
    }#i 
  }#if
}#t

### ====    6. GENERATE DEAD RECOVERIES : y.dead[i,t] & y.dead.binary[i,t] ====
## list recovered individuals
z.dead <- apply(true.z, 2, function(x){which(x == 3)})

## Initialize matrix of dead recoveries
y.dead <- matrix(0, nrow = dim(true.z)[1], ncol = dim(true.z)[2])

## Check that at least one individual was recovered
if(length(z.dead) > 0){
  ## Sample individual detections
  for(t in 1:n.years){
    #--- Check that at least one individual was recovered
    if(length(z.dead[[t]]) > 0){
      for(i in z.dead[[t]]){
        #--- Calculate squared distance to all detectors
        D2 <- (detector.xy[ ,1] - sCoords[i,1,t])^2 + (detector.xy[ ,2] - sCoords[i,2,t])^2
        #--- Calculate detector-specific (relative) recovery probability 
        P <- exp(-D2/(2 * sim$sigma * sim$sigma))
        #--- Sample individual recovery location 
        y.dead[i,t] <- rcat(n = 1, prob = P)
      }#i 
    }#if
  }
}#t

y.dead.binary <- y.dead
y.dead.binary[y.dead.binary > 0] <- 1

### ====    7. AUGMENT DATA ====
y.alive <- abind( y.alive, array(0, c(dim(y.alive)[1]/2,dim(y.alive)[2:3])), along = 1)
y.dead <- rbind(y.dead, matrix(0, nrow = nrow(y.dead)/2, ncol = ncol(y.dead)))
y.dead.binary <- rbind(y.dead.binary, matrix(0, nrow = nrow(y.dead.binary)/2, ncol = ncol(y.dead.binary)))
true.z.aug <- rbind(true.z, matrix(1, nrow = nrow(true.z)/2, ncol = ncol(true.z)))

## ---------------------------------------------------------------------------------
## ------ III.OPSCR MODEL FITTING ------
### ====    1.GENERATE RANDOM INITIAL VALUES ====
y.state <- apply(y.alive, c(1,3), function(x)ifelse(any(x > 0),2,1))

STATE <- matrix(c(1,1,0,
                  0,1,1,
                  0,0,1), nrow = 3, byrow = T)

OBS <- matrix(c(1,0,
                1,1,
                1,0), nrow = 3, byrow = T)

z.mx <- MakeZ( y = y.state,
               STATE = STATE,
               OBSERVATION = OBS,
               f.state = c(1,2),
               z.init.NA = TRUE)
z.mx$z.init.mx[(dim(true.z)[1]+1):dim(z.mx$z.init.mx)[1], ] <- 1

sxy.init <- MakeSxyInits( y = y.alive,
                          detCoords = as.matrix(detector.xy),
                          habitatPoly = bufferedStudyArea,
                          maxDetDist = 4*sim$sigma,
                          maxMoveDist = 4*sim$tau,
                          plot.check = F)

OPSCR_nimInits <- list( z = z.mx$z.init.mx,
                        sxy = sxy.init,
                        tau = runif(1,2,10),
                        sigma = runif(1,0,5),
                        p0 = runif(1,0,0.5),
                        phi = runif(1,0.2,0.8),
                        rho = runif(1,0.2,0.8),
                        gamma0 = runif(1,0,0.75))

### ====    2.NIMBLE PARAMETERS ====
OPSCR_nimParams <- c("N","gamma0","rho","p0","phi","sigma","tau")

### ====    3.NIMBLE CONSTANTS ====
OPSCR_nimConstants <- list( n.cells = numHabCells,
                            n.individuals = dim(y.alive)[1],
                            n.detectors = dim(y.alive)[2],
                            n.years = dim(y.alive)[3],
                            n.years1 = dim(y.alive)[3]-1)

### ====    4.NIMBLE DATA ====
OPSCR_nimData <- list( z = z.mx$z.reconstruct.mx,
                       y.alive = y.alive,
                       lowerHabCoords = lowerCoords,
                       upperHabCoords = upperCoords,
                       mu = rep(1, numHabCells),
                       detector.xy = detector.xy)

### ====    5.NIMBLE MODEL ====
OPSCR_nimModel <- nimbleCode({
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  tau ~ dgamma(0.001,0.001)
  
  for(i in 1:n.individuals){
    sxy[i,1:2,1] ~ dbinomPPSingle( lowerHabCoords[1:n.cells,1:2],
                                   upperHabCoords[1:n.cells,1:2],
                                   mu[1:n.cells],
                                   1, n.cells)
    for(t in 2:n.years){
      sxy[i,1:2,t] ~ dbinomMNormSourcePPSingle( lowerHabCoords[1:n.cells,1:2],
                                                upperHabCoords[1:n.cells,1:2],
                                                sxy[i,1:2,t-1],
                                                tau,
                                                mu[1:n.cells],
                                                1, n.cells, -1)
    }#t
  }#i
  
  ##-----------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##     
  gamma0 ~ dunif(0,1)                                                  
  phi ~ dunif(0,1)  
  rho ~ dunif(0,5)
  
  omeg1[1:3] <- c(1-gamma0, gamma0, 0)   
  
  for(t in 2:n.years){
    gamma[t-1] <- (N[t-1] * rho) / n.available[t-1]
  }#t
  
  for(t in 1:n.years1){
    omega[1,1:3,t] <- c(1-gamma[t], gamma[t], 0)
    omega[2,1:3,t] <- c(0, phi, 1-phi)
    omega[3,1:3,t] <- c(0,0,1)
  }#t
  
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:3]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
    }#t 										 
  }#i 
  
  ##-----------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------## 
  sigma ~ dunif(0,100)
  p0 ~ dunif(0,1)
  for(i in 1:n.individuals){
    for(t in 1:n.years){
      d2[i,1:n.detectors,t] <- calculateDistance( sxy[i,1:2,t],
                                                  detector.xy[1:n.detectors,1:2],
                                                  z[i,t] == 2)
      
      y.alive[i,1:n.detectors,t] ~ dbern_LESS( pZero = p0,
                                               sigma = sigma,
                                               d2 = d2[i,1:n.detectors,t],
                                               maxDist = 100,
                                               indicator = (z[i,t]==2))
    }#t
  }#i
  
  ##----------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(t in 1:n.years){
    N[t] <- sum(z[1:n.individuals, t]==2) 
    n.available[t] <- sum(z[1:n.individuals, t]==1)
  }#t
}) 

### ====    6.NIMBLE RUNS ====
OPSCR_model <- nimbleModel(code = OPSCR_nimModel,
                           constants = OPSCR_nimConstants,
                           inits = OPSCR_nimInits,
                           data = OPSCR_nimData,
                           check = FALSE)
OPSCR_cmodel <- compileNimble(OPSCR_model)
OPSCR_cmodel$calculate()
OPSCR_ModelMCMCConf <- configureMCMC(OPSCR_model, monitors = OPSCR_nimParams, thin = 1)
OPSCR_ModelMCMC <- buildMCMC(OPSCR_ModelMCMCConf)
OPSCR_ModelMCMCComp <- compileNimble(OPSCR_ModelMCMC, project = OPSCR_cmodel)

## TAKES APPROX. 2.5 HOURS TO RUN 2 CHAINS OF 40000 ITERATIONS ON A REGULAR LAPTOP
OPSCR_Runtime <- system.time(OPSCR_nimOutput <- runMCMC( OPSCR_ModelMCMCComp,
                                                         niter = 4000,
                                                         nburnin = 1000,
                                                         nchains = 2,
                                                         samplesAsCodaMCMC = TRUE,
                                                         summary = TRUE))

### ====    7.COMPARE ESTIMATES & SIMULATED VALUES ====
OPSCR_results <- as.data.frame(OPSCR_nimOutput$summary$all.chains)

OPSCR_simValues <- list( "N[1]" = sum(true.z[ ,1] == 2),
                            "N[2]" = sum(true.z[ ,2] == 2),
                            "N[3]" = sum(true.z[ ,3] == 2),
                            "N[4]" = sum(true.z[ ,4] == 2),
                            "N[5]" = sum(true.z[ ,5] == 2),
                            "gamma0" = sim$N0/dim(OPSCR_nimData$y.alive)[1],
                            "p0" = sim$p0,
                            "phi" = sim$phi,
                            "rho" = sim$rep,
                            "sigma" = sim$sigma,
                            "tau" = sim$tau)

OPSCR_results$RB <- (OPSCR_results$Mean - unlist(OPSCR_simValues))/unlist(OPSCR_simValues)
OPSCR_results$CV <- OPSCR_results$St.Dev./OPSCR_results$Mean


for(p in 1:length(OPSCR_simValues)){
  par(mfrow = c(1,2))
  traceplot(OPSCR_nimOutput$samples[ ,names(OPSCR_simValues)[p]])
  abline(h = OPSCR_simValues[p], col = "black", lwd = 3, lty = 2)
  plot(density(unlist(OPSCR_nimOutput$samples[ ,names(OPSCR_simValues)[p]])), main = names(OPSCR_simValues)[p])
  abline(v = OPSCR_simValues[p], col = "red", lwd = 2)
}#p

## ---------------------------------------------------------------------------------
## ------ IV.OPSCR + DR MODEL FITTING ------
### ====    1.GENERATE RANDOM INITIAL VALUES ====
y.state <- apply(y.alive, c(1,3), function(x)ifelse(any(x > 0),1,0))
y.state <- y.state + ifelse(y.dead > 0, 3, 1)

STATE <- matrix(c(1,1,0,0,
                  0,1,1,1,
                  0,0,0,1,
                  0,0,0,1), nrow = 4, byrow = T)

OBS <- matrix(c(1,0,0,
                1,1,0,
                0,0,1,
                1,0,0), nrow = 4, byrow = T)

z.mx <- MakeZ( y = y.state,
               STATE = STATE,
               OBSERVATION = OBS,
               f.state = c(1,2),
               z.init.NA = TRUE)
z.mx$z.init.mx[(dim(true.z)[1]+1):dim(z.mx$z.init.mx)[1], ] <- 1

sxy.init  <- MakeSxyInits( y = y.alive,
                           detCoords = as.matrix(detector.xy),
                           habitatPoly = bufferedStudyArea,
                           maxDetDist = 4*sim$sigma,
                           maxMoveDist = 4*sim$tau,
                           plot.check = F)

OPSCR_DR_nimInits <- list( z = z.mx$z.init.mx,
                           sxy = sxy.init,
                           tau = runif(1,2,10),
                           sigma = runif(1,0,5),
                           r = runif(1,0,1),
                           p0 = runif(1,0,0.5),
                           phi = runif(1,0.2,0.8),
                           rho = runif(1,0.2,0.8),
                           gamma0 = runif(1,0,0.75))

### ====    2.NIMBLE PARAMETERS ====
OPSCR_DR_nimParams <- c("N","gamma0","rho","p0","phi","sigma","tau","r")

### ====    3.NIMBLE CONSTANTS ====
OPSCR_DR_nimConstants <- list( n.cells = numHabCells,
                               n.individuals = dim(y.alive)[1],
                               n.detectors = dim(y.alive)[2],
                               n.years = dim(y.alive)[3],
                               n.years1 = dim(y.alive)[3]-1)

### ====    4.NIMBLE DATA ====
OPSCR_DR_nimData <- list( z = z.mx$z.reconstruct.mx,
                          y.alive = y.alive,
                          y.dead = y.dead.binary,
                          lowerHabCoords = lowerCoords,
                          upperHabCoords = upperCoords,
                          mu = rep(1,numHabCells),
                          detector.xy = detector.xy)

### ====    5.NIMBLE MODEL ====
OPSCR_DR_nimModel <- nimbleCode({
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  tau ~ dgamma(0.001,0.001)
  
  for(i in 1:n.individuals){
    sxy[i,1:2,1] ~ dbinomPPSingle( lowerHabCoords[1:n.cells,1:2],
                                   upperHabCoords[1:n.cells,1:2],
                                   mu[1:n.cells],
                                   1, n.cells)
    for(t in 2:n.years){
      sxy[i,1:2,t] ~ dbinomMNormSourcePPSingle( lowerHabCoords[1:n.cells,1:2],
                                                upperHabCoords[1:n.cells,1:2],
                                                sxy[i,1:2,t-1],
                                                tau,
                                                mu[1:n.cells],
                                                1, n.cells, -1)
    }#t
  }#i
  
  ##-----------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  gamma0 ~ dunif(0,1)                                                  
  phi ~ dunif(0,1)  
  rho ~ dunif(0,2)
  r ~ dunif(0,1)
  omeg1[1:3] <- c(1-gamma0, gamma0, 0)   
  
  for(t in 2:n.years){
    gamma[t-1] <- (N[t-1] * rho) / n.available[t-1]
  }#t
  
  for(t in 1:n.years1){
    omega[1,1:4,t] <- c(1-gamma[t],gamma[t],0,0)
    omega[2,1:4,t] <- c(0,phi,(1-phi)*r,(1-phi)*(1-r))
    omega[3,1:4,t] <- c(0,0,0,1)
    omega[4,1:4,t] <- c(0,0,0,1)
  }#t
  
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:3]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:4,t]) 
    }#t 										 
  }#i 
  
  ##-----------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------## 
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  for(i in 1:n.individuals){
    for(t in 1:n.years){
      d2[i,1:n.detectors,t] <- calculateDistance( sxy[i,1:2,t],
                                                  detector.xy[1:n.detectors,1:2])
      
      y.alive[i,1:n.detectors,t] ~ dbern_LESS( pZero = p0,
                                               sigma = sigma,
                                               d2 = d2[i,1:n.detectors,t],
                                               maxDist = 100,
                                               indicator = (z[i,t]==2))
      
      y.dead[i,t] ~ dbern(z[i,t]==3)
    }#t
  }#i
  ##----------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(t in 1:n.years){
    N[t] <- sum(z[1:n.individuals, t]==2) 
    n.available[t] <- sum(z[1:n.individuals,t]==1)
  }#t
})    

### ====    6.NIMBLE RUNS ====
OPSCR_DR_model <- nimbleModel(code = OPSCR_DR_nimModel,
                              constants = OPSCR_DR_nimConstants,
                              inits = OPSCR_DR_nimInits,
                              data = OPSCR_DR_nimData,
                              check = FALSE)
OPSCR_DR_cmodel <- compileNimble(OPSCR_DR_model)
OPSCR_DR_cmodel$calculate()
OPSCR_DR_ModelMCMCConf <- configureMCMC(OPSCR_DR_model, monitors = OPSCR_DR_nimParams, thin = 1)
OPSCR_DR_ModelMCMC <- buildMCMC(OPSCR_DR_ModelMCMCConf)
OPSCR_DR_ModelMCMCComp <- compileNimble(OPSCR_DR_ModelMCMC, project = OPSCR_DR_cmodel)

## TAKES APPROX. 2.5 HOURS TO RUN 2 CHAINS OF 40000 ITERATIONS ON A REGULAR LAPTOP
OPSCR_DR_Runtime <- system.time(OPSCR_DR_nimOutput <- runMCMC( OPSCR_DR_ModelMCMCComp,
                                                               niter = 4000,
                                                               nburnin = 1000,
                                                               nchains = 2,
                                                               samplesAsCodaMCMC = TRUE,
                                                               summary = TRUE))

### ====    7.COMPARE ESTIMATES & SIMULATED VALUES ====
OPSCR_DR_results <- as.data.frame(OPSCR_DR_nimOutput$summary$all.chains)

OPSCR_DR_simValues <- list( "N[1]" = sum(true.z[ ,1] == 2),
                          "N[2]" = sum(true.z[ ,2] == 2),
                          "N[3]" = sum(true.z[ ,3] == 2),
                          "N[4]" = sum(true.z[ ,4] == 2),
                          "N[5]" = sum(true.z[ ,5] == 2),
                          "gamma0" = sim$N0/dim(OPSCR_DR_nimData$y.alive)[1],
                          "p0" = sim$p0,
                          "phi" = sim$phi,
                          "r" = sim$r,
                          "rho" = sim$rep,
                          "sigma" = sim$sigma,
                          "tau" = sim$tau)

OPSCR_DR_results$RB <- (OPSCR_DR_results$Mean - unlist(OPSCR_DR_simValues))/unlist(OPSCR_DR_simValues)
OPSCR_DR_results$CV <- OPSCR_DR_results$St.Dev./OPSCR_DR_results$Mean

for(p in 1:length(OPSCR_DR_simValues)){
  par(mfrow = c(1,2))
  traceplot(OPSCR_DR_nimOutput$samples[ ,names(OPSCR_DR_simValues)[p]])
  abline(h = OPSCR_DR_simValues[p], col = "black", lwd = 3, lty = 2)
  plot(density(unlist(OPSCR_DR_nimOutput$samples[ ,names(OPSCR_DR_simValues)[p]])), main = names(OPSCR_DR_simValues)[p])
  abline(v = OPSCR_DR_simValues[p], col = "red", lwd = 2)
}#p

## ---------------------------------------------------------------------------------
## ------ V.OPSC2R MODEL FITTING ------
### ====    1.GENERATE RANDOM INITIAL VALUES ====
y.state <- apply(y.alive, c(1,3), function(x)ifelse(any(x > 0),1,0))
y.state <- y.state + ifelse(y.dead > 0, 3, 1)

STATE <- matrix(c(1,1,0,0,
                  0,1,1,1,
                  0,0,0,1,
                  0,0,0,1), nrow = 4, byrow = T)

OBS <- matrix(c(1,0,0,
                1,1,0,
                0,0,1,
                1,0,0), nrow = 4, byrow = T)

z.mx <- MakeZ( y = y.state,
               STATE = STATE,
               OBSERVATION = OBS,
               f.state = c(1,2),
               z.init.NA = TRUE)
z.mx$z.init.mx[(dim(true.z)[1]+1):dim(z.mx$z.init.mx)[1], ] <- 1

sxy.init  <- MakeSxyInits( y = y.alive,
                           detCoords = as.matrix(detector.xy),
                           habitatPoly = bufferedStudyArea,
                           maxDetDist = 4*sim$sigma,
                           maxMoveDist = 4*sim$tau,
                           plot.check = F)

OPSC2R_nimInits <- list( z = z.mx$z.init.mx,
                         sxy = sxy.init,
                         tau = runif(1,2,10),
                         sigma = runif(1,0,5),
                         r = runif(1,0,1),
                         p0 = runif(1,0,0.5),
                         phi = runif(1,0.2,0.8),
                         rho = runif(1,0.2,0.8),
                         gamma0 = runif(1,0,0.75))

### ====    2.NIMBLE PARAMETERS ====
OPSC2R_nimParams <- c("N","gamma0","rho","p0","phi","sigma","tau","r")

### ====    3.NIMBLE CONSTANTS ====
OPSC2R_nimConstants <- list( n.cells = numHabCells,
                             n.individuals = dim(y.alive)[1],
                             n.detectors = dim(y.alive)[2],
                             n.years = dim(y.alive)[3],
                             n.years1 = dim(y.alive)[3]-1)

### ====    4.NIMBLE DATA ====
OPSC2R_nimData <- list( z = z.mx$z.reconstruct.mx,
                        y.alive = y.alive,
                        y.dead = y.dead,
                        lowerHabCoords = lowerCoords,
                        upperHabCoords = upperCoords,
                        mu = rep(1,numHabCells),
                        detector.xy = detector.xy)

### ====    5.NIMBLE MODEL ====
OPSC2R_nimModel <- nimbleCode({
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  tau ~ dgamma(0.001,0.001)
  
  for(i in 1:n.individuals){
    sxy[i,1:2,1] ~ dbinomPPSingle( lowerHabCoords[1:n.cells,1:2],
                                   upperHabCoords[1:n.cells,1:2],
                                   mu[1:n.cells],
                                   1, n.cells)
    for(t in 2:n.years){
      sxy[i,1:2,t] ~ dbinomMNormSourcePPSingle( lowerHabCoords[1:n.cells,1:2],
                                                upperHabCoords[1:n.cells,1:2],
                                                sxy[i,1:2,t-1],
                                                tau,
                                                mu[1:n.cells],
                                                1, n.cells, -1)
    }#t
  }#i
  
  ##-----------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  gamma0 ~ dunif(0,1)                                                  
  phi ~ dunif(0,1)  
  rho ~ dunif(0,2)
  r ~ dunif(0,1)
  omeg1[1:3] <- c(1-gamma0, gamma0, 0)   
  
  ## Calculate gamma from rho as a function of n.available.
  for(t in 2:n.years){
    gamma[t-1] <- (N[t-1] * rho) / n.available[t-1]
  }#t
  
  for(t in 1:n.years1){
    omega[1,1:4,t] <- c(1-gamma[t],gamma[t],0,0)
    omega[2,1:4,t] <- c(0,phi,(1-phi)*r,(1-phi)*(1-r))
    omega[3,1:4,t] <- c(0,0,0,1)
    omega[4,1:4,t] <- c(0,0,0,1)
  }#t
  
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:3]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:4,t]) 
    }#t 										 
  }#i 
  
  ##-----------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------## 
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  for(i in 1:n.individuals){
    for(t in 1:n.years){
      d2[i,1:n.detectors,t] <- calculateDistance( sxy[i,1:2,t],
                                                  detector.xy[1:n.detectors,1:2])
      
      y.alive[i,1:n.detectors,t] ~ dbern_LESS( pZero = p0,
                                               sigma = sigma,
                                               d2 = d2[i,1:n.detectors,t],
                                               maxDist = 100,
                                               indicator = (z[i,t]==2))
      
      y.dead[i,t] ~ dcat_LESS( pZero =  1,
                               sigma = sigma,
                               d2 = d2[i,1:n.detectors,t],
                               maxDist = 100,
                               indicator = (z[i,t]==3))
    }#t
  }#i
  ##----------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(t in 1:n.years){
    N[t] <- sum(z[1:n.individuals, t]==2) 
    n.available[t] <- sum(z[1:n.individuals,t]==1)
  }#t
})   

### ====    6.NIMBLE RUNS ====
OPSC2R_model <- nimbleModel(code = OPSC2R_nimModel,
                            constants = OPSC2R_nimConstants,
                            inits = OPSC2R_nimInits,
                            data = OPSC2R_nimData,
                            check = FALSE)
OPSC2R_cmodel <- compileNimble(OPSC2R_model)
OPSC2R_cmodel$calculate()
OPSC2R_ModelMCMCConf <- configureMCMC(OPSC2R_model, monitors = OPSC2R_nimParams, thin = 1)
OPSC2R_ModelMCMC <- buildMCMC(OPSC2R_ModelMCMCConf)
OPSC2R_ModelMCMCComp <- compileNimble(OPSC2R_ModelMCMC, project = OPSC2R_cmodel)

## TAKES APPROX. 2.5 HOURS TO RUN 2 CHAINS OF 40000 ITERATIONS ON A REGULAR LAPTOP
OPSC2R_Runtime <- system.time(OPSC2R_nimOutput <- runMCMC( OPSC2R_ModelMCMCComp,
                                                    niter = 4000,
                                                    nburnin = 1000,
                                                    nchains = 2,
                                                    samplesAsCodaMCMC = TRUE,
                                                    summary = TRUE))

### ====    7.COMPARE ESTIMATES & SIMULATED VALUES ====
OPSC2R_results <- as.data.frame(OPSC2R_nimOutput$summary$all.chains)

OPSC2R_simValues <- list( "N[1]" = sum(true.z[ ,1] == 2),
                "N[2]" = sum(true.z[ ,2] == 2),
                "N[3]" = sum(true.z[ ,3] == 2),
                "N[4]" = sum(true.z[ ,4] == 2),
                "N[5]" = sum(true.z[ ,5] == 2),
                "gamma0" = sim$N0/dim(OPSC2R_nimData$y.alive)[1],
                "p0" = sim$p0,
                "phi" = sim$phi,
                "r" = sim$r,
                "rho" = sim$rep,
                "sigma" = sim$sigma,
                "tau" = sim$tau)
                
OPSC2R_results$RB <- (OPSC2R_results$Mean - unlist(OPSC2R_simValues))/unlist(OPSC2R_simValues)
OPSC2R_results$CV <- OPSC2R_results$St.Dev./OPSC2R_results$Mean

for(p in 1:length(OPSC2R_simValues)){
  par(mfrow = c(1,2))
  traceplot(OPSC2R_nimOutput$samples[ ,names(OPSC2R_simValues)[p]])
  abline(h = OPSC2R_simValues[p], col = "black", lwd = 3, lty = 2)
  plot(density(unlist(OPSC2R_nimOutput$samples[ ,names(OPSC2R_simValues)[p]])), main = names(OPSC2R_simValues)[p])
  abline(v = OPSC2R_simValues[p], col = "red", lwd = 2)
}#p

## ---------------------------------------------------------------------------------
## ------ VI.COMPARE MODEL RESULTS ------
OPSCR_results
OPSCR_DR_results
OPSC2R_results
## ---------------------------------------------------------------------------------
