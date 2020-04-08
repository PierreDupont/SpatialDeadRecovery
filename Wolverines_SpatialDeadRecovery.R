#############################################################################
####### ---- Open-Population Spatial Capture-Recapture-Recovery ---- ########
####### ------------------ WOLVERINE DATA ANALYSIS ----------------- ########
#############################################################################
## ------ IMPORT REQUIRED LIBRARIES ------
rm(list=ls())
library(nimble)

## ------ SOURCE THE REQUIRED FUNCTIONS ------
source("FUNCTIONS/dbin_LESS_Cached_MultipleCov.R")
source("FUNCTIONS/dcat_Cached_Sparse.R")
source("FUNCTIONS/pointProcess.R")

## ----------------------------------------------------------------------------------------------
## ------ I.LOAD WOLVERINE DATASETS -----
## FEMALE DATASETS
load("DATA/female_OPSCR_data.RData")
load("DATA/female_OPSC2R_data.RData")

## MALE DATASETS
# load("DATA/male_OPSCR_data.RData")
# load("DATA/male_OPSC2R_data.RData")

## ----------------------------------------------------------------------------------------------
## ------ II.MODELS DEFINITION ------- 
## OPSC2R model
wolverine_OPSC2R_model <- nimbleCode({
  ##--------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  tau ~ dunif(0,10)                                                      # movement parameter
  betaDens  ~ dnorm(0.0,0.01)                                            # Cooefficient of the number of dens on density
    mu[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])  # Expected density per habitat cell
  
  for(i in 1:n.individuals){
    # Initial AC placement (= spatial point process with intensity mu[1:numHabWindows])
    sxy[i,1:2,1] ~ dbinomPPSingle( lowerHabCoords[1:numHabWindows,1:2],
                                   upperHabCoords[1:numHabWindows,1:2],
                                   mu[1:numHabWindows],
                                   1, numHabWindows)
    
    for(t in 2:n.years){
      # AC movement(= bivariate normal point process centered on sxy[i,1:2,t-1]
      # weighted by intensity mu[1:numHabWindows])
      sxy[i,1:2,t] ~ dbinomMNormSourcePPSingle( lowerHabCoords[1:numHabWindows,1:2],
                                                upperHabCoords[1:numHabWindows,1:2],
                                                sxy[i,1:2,t-1],
                                                tau,
                                                mu[1:numHabWindows],
                                                1, numHabWindows, -1)
    }#t
  }#i
  
  ##--------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##
  omeg1[1:2] ~ ddirch(alpha[1:2]) # Initial state assignment ("available" or "alive")  
  
  for(t in 1:n.years1){
    # PRIORS 
    gamma[t] ~ dunif(0,1)         # recruitment probability
    w[t] ~ dunif(0,0.59)          # "natural" mortality (= (1-phi) * (1-r) )
    h[t] ~ dunif(0,0.39)          # human-caused mortality (= (1-phi) * r )
    phi[t] <- 1-h[t]-w[t]         # survival probability 
    
    # state transitions for individuals in state "AVAILABLE"
    omega[1,1:4,t] <- c(1-gamma[t], gamma[t], 0, 0)
    # state transitions for individuals in state "ALIVE"
    omega[2,1:4,t] <- c(0, phi[t], h[t], w[t])     
    # state transitions for individuals in state "RECOVERED"
    omega[3,1:4,t] <- c(0, 0, 0, 1)
    # state transitions for individuals in state "DEAD"
    omega[4,1:4,t] <- c(0, 0, 0, 1)
  }#t
  
  ## STATE TRANSITIONS
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:2]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:4,t]) 
    }#i 								
  }#t 
  
    ##---------------------------------------------------------------------------------------------   
    ##-----------------------------##
    ##----- DETECTION PROCESS -----## 
    ##-----------------------------##
    sigma ~ dunif(0,4)          # scale parameter of the detection function
    betaResponse ~ dunif(-5,5)  # effect of previous detection on p0
    
    for(c in 1:n.covs){         # betaCovs[1] : effect of GPS tracks on p0
      betaCovs[c] ~ dunif(-5,5) # betaCovs[2] : effect of distance to roads on p0
    }#c                         # betaCovs[3] : effect of snow cover on p0
    
    for(c in 1:n.counties){
      for(t in 1:n.years){
        p0[c,t] ~ dunif(0,1)    # County and year-specific baseline detection probabilities
      }#t
    }#c  
    
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){
      # ALIVE DETECTIONS
      y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCov( sxy = sxy[i,1:2,t],
                                                                   sigma = sigma,
                                                                   p0 = p0[1:n.counties,t],
                                                                   nbDetections = nbDetections[i,t],
                                                                   yDets = yDets[i,1:nMaxDetectors,t],
                                                                   detector.xy = detector.xy[1:n.detectors,1:2],
                                                                   trials = trials[1:n.detectors],
                                                                   detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                   nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                   resizeFactor = resizeFactor,
                                                                   maxNBDets = maxNBDets,
                                                                   habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                   indicator = (z[i,t]==2),
                                                                   detCounties = detCounties[1:n.detectors],
                                                                   detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                   betaCov = betaCovs[1:n.covs],
                                                                   betaResponse = betaResponse,
                                                                   detResponse = detResponse[i,t])
      # DEAD RECOVERY
      y.dead[i,t] ~ dcat_Cached_Sparse( pZero = 1,
                                        sigma = sigma,
                                        sxy = sxy[i,1:2,t],
                                        detectorCoords = detector.dead.xy[1:n.detectors.dead,1:2],
                                        detectorID = detectorIndex.dead[1:n.cellsSparse.dead,1:maxNBDets.dead],
                                        detectorNum = nDetectorsLESS.dead[1:n.cellsSparse.dead],
                                        habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet], 
                                        habitatFactor = resizeFactor,
                                        indicator = (z[i,t] == 3))
    }#i
  }#t
  
  ##---------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(t in 1:n.years){
    N[t] <- sum(z[1:n.individuals,t] == 2) # POPULATION SIZE
  }#t
})

## OPSCR model
wolverine_OPSCR_model <- nimbleCode({
  ##--------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  tau ~ dunif(0,10)                                                      # movement parameter
  betaDens  ~ dnorm(0.0,0.01)                                            # Cooefficient of the number of dens on density
  mu[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])  # Expected density per habitat cell
  
  for(i in 1:n.individuals){
    # Initial AC placement (= spatial point process with intensity mu[1:numHabWindows])
    sxy[i,1:2,1] ~ dbinomPPSingle( lowerHabCoords[1:numHabWindows,1:2],
                                   upperHabCoords[1:numHabWindows,1:2],
                                   mu[1:numHabWindows],
                                   1, numHabWindows)
    
    for(t in 2:n.years){
      # AC movement(= bivariate normal point process centered on sxy[i,1:2,t-1]
      # weighted by intensity mu[1:numHabWindows])
      sxy[i,1:2,t] ~ dbinomMNormSourcePPSingle( lowerHabCoords[1:numHabWindows,1:2],
                                                upperHabCoords[1:numHabWindows,1:2],
                                                sxy[i,1:2,t-1],
                                                tau,
                                                mu[1:numHabWindows],
                                                1, numHabWindows, -1)
    }#t
  }#i
  
  ##--------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##
  omeg1[1:2] ~ ddirch(alpha[1:2]) # Initial state assignment ("available" or "alive")  
  
  for(t in 1:n.years1){
    # PRIORS 
    gamma[t] ~ dunif(0,1)         # recruitment probability
    w[t] ~ dunif(0,1)             # "natural" mortality (= (1-phi) * (1-r) )
    phi[t] <- 1-w[t]              # survival probability 
    
    # state transitions for individuals in state "AVAILABLE"
    omega[1,1:3,t] <- c(1-gamma[t], gamma[t], 0)
    # state transitions for individuals in state "ALIVE"
    omega[2,1:3,t] <- c(0, phi[t], w[t])     
    # state transitions for individuals in state "DEAD"
    omega[3,1:3,t] <- c(0, 0, 1)
  }#t
  
  ## STATE TRANSITIONS
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:2]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
    }#i 								
  }#t 
  
  ##---------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  sigma ~ dunif(0,4)          # scale parameter of the detection function
  betaResponse ~ dunif(-5,5)  # effect of previous detection on p0
  
  for(c in 1:n.covs){         # betaCovs[1] : effect of GPS tracks on p0
    betaCovs[c] ~ dunif(-5,5) # betaCovs[2] : effect of distance to roads on p0
  }#c                         # betaCovs[3] : effect of snow cover on p0

  for(c in 1:n.counties){
    for(t in 1:n.years){
      p0[c,t] ~ dunif(0,1)    # County and year-specific baseline detection probabilities
    }#t
  }#c  
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){
      # ALIVE DETECTIONS
      y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCov( sxy = sxy[i,1:2,t],
                                                                   sigma = sigma,
                                                                   p0 = p0[1:n.counties,t],
                                                                   nbDetections = nbDetections[i,t],
                                                                   yDets = yDets[i,1:nMaxDetectors,t],
                                                                   detector.xy = detector.xy[1:n.detectors,1:2],
                                                                   trials = trials[1:n.detectors],
                                                                   detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                   nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                   resizeFactor = resizeFactor,
                                                                   maxNBDets = maxNBDets,
                                                                   habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                   indicator = (z[i,t]==2),
                                                                   detCounties = detCounties[1:n.detectors],
                                                                   detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                   betaCov = betaCovs[1:n.covs],
                                                                   betaResponse = betaResponse,
                                                                   detResponse = detResponse[i,t])
    }#i
  }#t
  
  ##---------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(t in 1:n.years){
    N[t] <- sum(z[1:n.individuals,t] == 2)   # POPULATION SIZE
  }#t
})

## ----------------------------------------------------------------------------------------------
## ------ III.NIMBLE MODELS FITTING -----
## RUN OPSC2R model
OPSC2R_model <- nimbleModel(code = wolverine_OPSC2R_model,
                           constants = wolverine_OPSC2R_constants,
                           inits = wolverine_OPSC2R_inits,
                           data = wolverine_OPSC2R_data,
                           check = FALSE, calculate = FALSE)
OPSC2R_cmodel <- compileNimble(OPSC2R_model)
OPSC2R_cmodel$calculate()
OPSC2R_ModelMCMCConf <- configureMCMC(OPSC2R_model, monitors = wolverine_OPSC2R_params, thin = 1)
OPSC2R_ModelMCMC <- buildMCMC(OPSC2R_ModelMCMCConf)
OPSC2R_ModelMCMCComp <- compileNimble(OPSC2R_ModelMCMC, project = OPSC2R_cmodel)

## TAKES APPROX.  HOURS TO RUN 1000 ITERATIONS ON A REGULAR LAPTOP
OPSC2R_Runtime <- system.time(wolverine_OPSC2R_results <- runMCMC( OPSC2R_ModelMCMCComp,
                                                                   niter = 1000,
                                                                   nburnin = 0,
                                                                   nchains = 1,
                                                                   samplesAsCodaMCMC = TRUE))


## RUN OPSCR model
OPSCR_model <- nimbleModel(code = wolverine_OPSCR_model,
                           constants = wolverine_OPSCR_constants,
                           inits = wolverine_OPSCR_inits,
                           data = wolverine_OPSCR_data,
                           check = FALSE, calculate = FALSE)
OPSCR_cmodel <- compileNimble(OPSCR_model)
OPSCR_cmodel$calculate()
OPSCR_ModelMCMCConf <- configureMCMC(OPSCR_model, monitors = wolverine_OPSCR_params, thin = 1)
OPSCR_ModelMCMC <- buildMCMC(OPSCR_ModelMCMCConf)
OPSCR_ModelMCMCComp <- compileNimble(OPSCR_ModelMCMC, project = OPSCR_cmodel)

## TAKES APPROX. 45 min TO RUN 1000 ITERATIONS ON A REGULAR LAPTOP
OPSCR_Runtime <- system.time(OPSCR_results <- runMCMC( OPSCR_ModelMCMCComp,
                                                       niter = 1000,
                                                       nburnin = 0,
                                                       nchains = 1,
                                                       samplesAsCodaMCMC = TRUE,
                                                       summary = TRUE))

## ----------------------------------------------------------------------------------------------
## ------ IV.LOAD PROCESSED WOLVERINE RESULTS ---- 
load("DATA/female_OPSCR_results.RData")
female_OPSCR_results$mean   ## posterior mean estimates (female_OPSCR_results$sims.list contains all MCMC samples)
female_OPSCR_results$q2.5   ## lower 95% CI 
female_OPSCR_results$q97.5  ## upper 95% CI
female_OPSCR_results$Rhat   ## Rhat values

load("DATA/female_OPSC2R_results.RData")
female_OPSC2R_results$mean
female_OPSC2R_results$q2.5
female_OPSC2R_results$q97.5
female_OPSC2R_results$Rhat

load("DATA/male_OPSCR_results.RData")
male_OPSCR_results$mean
male_OPSCR_results$q2.5
male_OPSCR_results$q97.5
male_OPSCR_results$Rhat

load("DATA/male_OPSC2R_results.RData")
male_OPSC2R_results$mean
male_OPSC2R_results$q2.5
male_OPSC2R_results$q97.5
male_OPSC2R_results$Rhat


## ----------------------------------------------------------------------------------------------