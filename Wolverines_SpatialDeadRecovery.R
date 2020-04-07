###################################################
##### ------------- WOLF : OPSCR ------------ #####
##### ------- AC[dbinomPP().dMNorm()] ------- #####
##### --- Z[gamma(t).phi(state.t).psi(t)] --- #####
##### ----- Y[p0(state.t).sigma(state)] ----- #####
###################################################
## ------ IMPORT REQUIRED LIBRARIES ------
rm(list=ls())
gc()
#.libPaths(new=c(.libPaths(),"C://PROJECTS//R//LIBRARY"))
library(rgdal)
library(raster)
library(sp)
library(coda)
library(nimble)
library(spdep)
library(rgeos)
library(maptools)
library(stringr)
library(abind)
library(R.utils)
library(adehabitatHR)
library(sf)
library(fasterize)
#library(snow)
            

## ------ SOURCE THE REQUIRED FUNCTIONS ------
sourceDirectory(dir.function, modifiedOnly = FALSE)
sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)
source("FUNCTIONS/dbin_LESS_Cached_MultipleCovResponse.R")
source("FUNCTIONS/dcat_Cached_Sparse.R")

## ----------------------------------------------------------------------------------------------
## ------ I.LOAD DATA -----
female_data_OPSCR <- load("DATA/female_data_OPSCR.RData")
# female_data_OPSC2R <- load("DATA/female_data_OPSC2R.RData")
# male_data_OPSCR <- load("DATA/male_data_OPSCR.RData")
# male_data_OPSC2R <- load("DATA/male_data_OPSC2R.RData")

## ----------------------------------------------------------------------------------------------
## ------ II.MODELS DEFINITION ------- 
wolverine_OPSC2R <- nimbleCode({
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
                                                                           nbDetections[i,t],
                                                                           yDets = yDets[i,1:nMaxDetectors,t],
                                                                           detector.xy = detector.xy[1:n.detectors,1:2],
                                                                           trials = trials[1:n.detectors],
                                                                           detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                           nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                           ResizeFactor = ResizeFactor,
                                                                           maxNBDets = maxNBDets,
                                                                           habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                           indicator = (z[i,t]==2),
                                                                           p0[1:n.countries,t],
                                                                           detCountries[1:n.detectors],
                                                                           detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                           betaCov = betaCovs[1:n.covs],
                                                                           BetaResponse = betaResponse,
                                                                           detResponse = detResponse[i,t])
      # DEAD RECOVERY
      y.dead[i,t] ~ dcat_Cached_Sparse( pZero = 1,
                                        sigma = sigma,
                                        sxy = sxy[i,1:2,t],
                                        detectorCoords = detector.dead.xy[1:n.detectors.dead,1:2],
                                        detectorID = detectorIndex.dead[1:n.cellsSparse.dead,1:maxNBDets.dead],
                                        detectorNum = nDetectorsLESS.dead[1:n.cellsSparse.dead],
                                        habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet], 
                                        habitatFactor = ResizeFactor,
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

wolverine_OPSCR <- nimbleCode({
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
    phi[t] <- 1-h[t]-w[t]         # survival probability 
    
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
      y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_Multiple( sxy = sxy[i,1:2,t],
                                                                sigma = sigma,
                                                                nbDetections[i,t],
                                                                yDets = yDets[i,1:nMaxDetectors,t],
                                                                detector.xy = detector.xy[1:n.detectors,1:2],
                                                                trials = trials[1:n.detectors],
                                                                detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                ResizeFactor = ResizeFactor,
                                                                maxNBDets = maxNBDets,
                                                                habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                indicator = (z[i,t]==2),
                                                                p0[1:n.counties,t],
                                                                detCountries[1:n.detectors],
                                                                detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                betaCov = betaCovs[1:n.covs],
                                                                BetaResponse = betaResponse,
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
## ------ III.OPSC2R MODEL FITTING
OPSC2R_model <- nimbleModel( code = wolverine_OPSC2R,
                      constants = wolverine_OPSC2R_constants,
                      inits = wolverine_OPSC2R_inits,
                      data = wolverine_OPSC2R_data,
                      check = FALSE,
                      calculate = FALSE) 
OPSC2R_model$calculate()  
OPSC2R_cmodel <- compileNimble(OPSC2R_model)
OPSC2R_conf <- configureMCMC(OPSC2R_model, monitors = wolverine_OPSC2R_params, thin = 1)
OPSC2R_Rmcmc <- buildMCMC(OPSC2R_conf)
OPSC2R_compiledList <- compileNimble(list(model = OPSC2R_model, mcmc = OPSC2R_Rmcmc))
OPSC2R_Cmodel <- OPSC2R_compiledList$model
OPSC2R_Cmcmc <- OPSC2R_compiledList$mcmc

## RUN NIMBLE MCMC IN SUCCESSIVE BITES DECREASE MEMORY USAGE
## SET NUMBER OF BITES AND NUMBER OF ITERATIONS PER BITE
bite.size <- 100 
bite.number <- 10

## LOOP OVER NUMBER OF BITES
for(nb in 1:bite.number){
  print(nb)
  if(nb == 1){
    ## run initial MCMC
    MCMCRuntime <- system.time(Cmcmc$run(bite.size))
  } else {      
    ## run subsequent MCMCs
    MCMCRuntime <- system.time(Cmcmc$run(bite.size, reset = FALSE))
  }
  
  ## STORE BITE OUTPUT IN A MATRIX
  mcmcSamples <- as.matrix(Cmcmc$mvSamples)
  CumulRunTime <- proc.time() - ptm
  
  ## EXPORT NIMBLE OUTPUT 
  outname <- file.path(path.OUT, paste("NimbleBite", nb, "_FOR", set, sep = ""))
  save(CumulRuntime, MCMCRuntime, mcmcSamples, file = outname)
  
  ## FREE UP MEMORY SPACE 
  rm("mcmcSamples") 
  Cmcmc$mvSamples$resize(0) ## reduce the internal mvSamples object to 0 rows,
  gc() ## run R's garbage collector
}#nb

OPSC22R_model <- nimbleModel(code = wolverine_OPSC2R_model,
                           constants = wolverine_OPSC2R_constants,
                           inits = wolverine_OPSC2R_inits,
                           data = wolverine_OPSC2R_data,
                           check = FALSE)
OPSC2R_cmodel <- compileNimble(OPSC2R_model)
OPSC2R_cmodel$calculate()
OPSC2R_ModelMCMCConf <- configureMCMC(OPSC2R_model, monitors = wolverine_OPSC2R_params, thin = 1)
OPSC2R_ModelMCMC <- buildMCMC(OPS2CR_ModelMCMCConf)
OPSC2R_ModelMCMCComp <- compileNimble(OPSC2R_ModelMCMC, project = OPSC2R_cmodel)

## TAKES APPROX. XXX HOURS TO RUN 5000 ITERATIONS ON A REGULAR LAPTOP
OPSC2R_Runtime <- system.time(wolverine_OPSC2R_output <- runMCMC( OPSC2R_ModelMCMCComp,
                                                         niter = 1000,
                                                         nburnin = 0,
                                                         nchains = 2,
                                                         samplesAsCodaMCMC = TRUE))

OPSCR_model <- nimbleModel(code = wolverine_OPSCR_model,
                           constants = wolverine_OPSCR_constants,
                           inits = wolverine_OPSCR_inits,
                           data = wolverine_OPSCR_data,
                           check = FALSE)
OPSCR_cmodel <- compileNimble(OPSCR_model)
OPSCR_cmodel$calculate()
OPSCR_ModelMCMCConf <- configureMCMC(OPSCR_model, monitors = OPSCR_nimParams, thin = 1)
OPSCR_ModelMCMC <- buildMCMC(OPSCR_ModelMCMCConf)
OPSCR_ModelMCMCComp <- compileNimble(OPSCR_ModelMCMC, project = OPSCR_cmodel)

## TAKES APPROX. XXX HOURS TO RUN 5000 ITERATIONS ON A REGULAR LAPTOP
OPSCR_Runtime <- system.time(OPSCR_nimOutput <- runMCMC( OPSCR_ModelMCMCComp,
                                                         niter = 1000,
                                                         nburnin = 0,
                                                         nchains = 2,
                                                         samplesAsCodaMCMC = TRUE,
                                                         summary = TRUE))
## ------ IV.PROCESS RESULTS ------
## -----------------------------------------------------------------------------------------------
### ==== 1. LOAD & PROCESS NIMBLE OUTPUTS ====
# List the directories containing bite outputs
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut", list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x)length(list.files(x))))
minBites <- min(numBites)

nimOutput <- list()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out <- list()
  for(x in 11:minBites){
    print(x)
    load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
    params.simple <- sapply(strsplit(colnames(mcmcSamples), "\\["), "[", 1)
    parmIndex <- which(! params.simple %in% c("sxy","z"))
    #thinIndex <- seq(1,dim(mcmcSamples)[1], 5)
    out[[x]] <- mcmcSamples[ ,parmIndex] 
  }#x
  out.mx <- do.call(rbind, out)
  nimOutput[[p]] <- as.mcmc(out.mx)
}#p

nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(nimOutput)#, params.omit = c("sxy","z"))
save(myResults, file = file.path(myVars$WD, myVars$modelName,"results_SDR.RData"))

## Process sxy & z posteriors
nimOutput_sxy <- list()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out <- list()
  for(x in 11:minBites){
    print(x)
    load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
    params.simple <- sapply(strsplit(colnames(mcmcSamples), "\\["), "[", 1)
    parmIndex <- which(params.simple %in% c("sxy","z"))
    thinIndex <- seq(1,dim(mcmcSamples)[1], 3)
    out[[x]] <- mcmcSamples[thinIndex,parmIndex] 
  }#x
  out.mx <- do.call(rbind, out)
  nimOutput_sxy[[p]] <- as.mcmc(out.mx)
}#p
nimOutput_sxy <- as.mcmc.list(nimOutput_sxy)
myResults_sxy <- ProcessCodaOutput(x = nimOutput_sxy)

## RESCALE THE COORDINATES 
sxyScaled <- GridToUTM(data.sxy = myResults_sxy$sims.list$sxy,
                       grid.sp = myHabitat$habitat.sp )
# rm(nimOutput)

### ==== 2. PLOT PARAMETERS ESTIMATES ====
### ====    2.1.N ==== 
### ====       2.1.1.By Country ====
habbR <- myHabitat$habitat.r

## Identify SWEDEN & NORWAY in the raster
SWE <- COUNTRIES[which(COUNTRIES$ISO %in% c("SWE")), ]     ## Just take Sweden
NOR <- COUNTRIES[which(COUNTRIES$ISO %in% c("NOR")), ]     ## Just take Norway
NOR <- gSimplify(spgeom = NOR, tol = 500)
SWE <- gSimplify(spgeom = SWE, tol = 500)
this.r <- RasterizePolygon( poly = SWE,
                            r = habbR,
                            CoverToKeepHabitat = 50,
                            fasterize = TRUE)
habbR[this.r == 1] <- 2
this.r <- RasterizePolygon( poly = NOR,
                            r = habbR,
                            CoverToKeepHabitat = 50,
                            fasterize = TRUE)
habbR[this.r == 1] <- 1
plot(habbR)

## Remove buffer from habitat
e <- extent(myStudyArea)
e.sp <- as(e, 'SpatialPolygons')  
habbR <- mask(habbR, e.sp)

## Convert to text 
habbR[habbR == 2] <- "SWE"
habbR[habbR == 1] <- "NOR"
habbR[habbR == 0] <- NA
gc()

## EXTRACT POPULATION SIZE per COUNTRY 
DensiCountries <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
  DensiCountries[[t]] <- EstimateN_v3( habRaster = habbR,
                                       posterior.sxy = sxyScaled$data.scaled.xy[,,,t],
                                       posterior.z = myResults_sxy$sims.list$z[,,t],
                                       alive.states = 2,
                                       return.all = FALSE,
                                       regionEstimates = TRUE)
  print(t)
}#t
DensiCountries[[1]]$PosteriorsRegion$`1`
DensiCountries[[1]]$summary

### ====       2.1.2.By County ====
habbRCounties <- myHabitat$habitat.r
plot(myHabitat$habitat.r)

## Get the counties 
COUNTIES <- aggregate(x = COMMUNES, by = "NAME_1")
COUNTIESsimp <- gSimplify(COUNTIES,tol = 500, topologyPreserve = T)
COUNTIESsimp$NAME_1 <- COUNTIES$NAME_1
plot(myHabitat$habitat.poly,add=T)
plot(COUNTIESsimp, add = T)

##remove buffer from habitat
habbRCountieswbuff <- habbRCounties <- mask(habbRCounties, e.sp)

for(i in 1:length(COUNTIESsimp)){
  # Identify SWEDEN in the raster
  this.r <- RasterizePolygon( poly = COUNTIESsimp[i,],
                              r = habbRCountieswbuff,
                              CoverToKeepHabitat = 50,
                              fasterize = TRUE)
  habbRCounties[this.r==1] <- as.character(COUNTIESsimp$NAME_1[i])
}
plot(habbRCounties)
plot(myHabitat$habitat.poly,add=T)

## EXTRACT POPULATION SIZE per COUNTY 
DensiCounties <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
  DensiCounties[[t]] <- EstimateN_v3( habRaster = habbRCounties, 
                                      posterior.sxy = sxyScaled$data.scaled.xy[,,,t],
                                      posterior.z = myResults_sxy$sims.list$z[,,t],
                                      alive.states = 2,
                                      return.all = FALSE,
                                      regionEstimates = TRUE)
  print(t)
}#t
## They should be the same. 
DensiCounties[[t]]$summary["Total", ]
DensiCountries[[t]]$summary["Total", ]

### ====       2.1.3.PLOTS ====
pdf(file = file.path( myVars$WD, myVars$modelName,
                      paste(myVars$modelName, "_NRegions.pdf", sep = "")))
## TOTAL 
plot(-1000, xlim = c(0.5,nYears+0.5), ylim = c(0,1000), ylab = "N Total", xaxt = "n")
axis(1, at = c(1:nYears), labels = years)
country.colors <- c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
for(t in 1:nYears){
  plot.violins2( DensiCountries[[t]]$PosteriorsRegion["ALL"],
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "red",
                 alpha = 0.2,
                 border.col = "red",
                 add = T)
  text( DensiCountries[[t]]$summary["Total", "median"],
        x = t, 
        y = DensiCountries[[t]]$summary["Total","median"] + 150)
}#t

## COUNTRIES 
plot(-1000, xlim = c(0.5, nYears + 0.5), ylim = c(0,700), ylab = "N", xaxt = "n")
axis(1, at = c(1:nYears), labels = years)
country.colors <- c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
for(t in 1:nYears){
  plot.violins2(DensiCountries[[t]]$PosteriorsRegion["NOR"],
                x = t,
                at = t - 0.15,
                violin.width = 0.3,
                col = country.colors[1],
                alpha = 0.2,
                border.col = country.colors[1],
                add = T)
  text( DensiCountries[[t]]$summary["NOR","median"],
        x = t - 0.15,
        y = DensiCountries[[t]]$summary["NOR","median"] - 50)
  
  plot.violins2(DensiCountries[[t]]$PosteriorsRegion["SWE"],
                x = t,
                at = t + 0.15,
                violin.width = 0.3,
                col = country.colors[2],
                alpha = 0.2,
                border.col =  country.colors[2],
                add = T)
  text( DensiCountries[[t]]$summary["SWE","median"],
        x = t + 0.15,
        y = DensiCountries[[t]]$summary["SWE","median"] + 50)
}#t
legend("topright", legend = names(country.colors), fill = country.colors, bty = "n")


## COUNTIES
idcounty <- rownames(DensiCounties[[1]]$summary)
par(mfrow = c(2,2))

for(co in 2:(length(idcounty)-1)){
  plot(-1000, xlim=c(0.5,nYears+0.5), ylim=c(0,300), ylab = paste("N", idcounty[co],sep=""),xaxt="n")
  axis(1, at = c(1:nYears), labels = years)
  for(t in 1:nYears){
    plot.violins2( DensiCounties[[t]]$PosteriorsRegion[idcounty[co]],
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = "firebrick4",
                   alpha = 0.2,
                   border.col = "firebrick4",
                   add = T)
    text( DensiCounties[[t]]$summary[idcounty[co],"median"],
          x = t,
          y = DensiCounties[[t]]$summary[idcounty[co],"median"] + 100)
  }#t      
  
  plot(myStudyArea)
  plot(COUNTIES, add = T)
  plot(COUNTIES[COUNTIES$NAME_1 == idcounty[co],], col = "red", add = T)
}#co
dev.off()      

### ====    2.2.ESTIMATES ====
pdf(file = file.path( myVars$WD, myVars$modelName,
                      paste(myVars$modelName, "_ESTIMATES.pdf", sep = "")))

### ====       2.2.1.N ====
par(mfrow = c(1,1), mar = c(5,5,5,5))
plot(10, xlim = c(0, nYears+1), ylim = c(200,800), type ="n", xaxt = "n", xlab = "Years", ylab = "N")
axis(1, c(1:nYears),labels = years)
for(t in 1:nYears){
  plot.violins(list(myResults$sims.list$N[,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = "firebrick3",
               add = T,
               alpha = 0.3,
               border.col = "firebrick3")
  text( x = t,
        y = mean(myResults$sims.list$N[,t]) + 100,
        labels = round(mean(myResults$sims.list$N[ ,t])))
}#t

params <- dimnames(nimOutput[[1]])[[2]][grep("N",dimnames(nimOutput[[1]])[[2]])]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}#i

### ====       2.2.2.h ==== 
par(mfrow = c(1,1))
plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "h")
axis(1, c(1:nYears),labels = years)
myCol <- "firebrick3"
for(t in 1:(nYears-1)){
  plot.violins(list(myResults$sims.list$h[ ,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t

params <- dimnames(nimOutput[[1]])[[2]][grep("h\\[",dimnames(nimOutput[[1]])[[2]])]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}

### ====       2.2.3.gamma ==== 
par(mfrow = c(1,1))
plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "gamma")
axis(1, c(1:nYears),labels = years)
myCol <- "firebrick3"
for(t in 1:(nYears-1)){
  plot.violins(list(myResults$sims.list$gamma[ ,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t

params <- dimnames(nimOutput[[1]])[[2]][grep("gamma",dimnames(nimOutput[[1]])[[2]])]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}

### ====       2.2.4.rho ==== 
posterior.rho <- t(apply(X = myResults_sxy$sims.list$z, MARGIN = 1, FUN = function(x){
  rho <- NULL
  for(t in 2:(dim(x)[2])){
    Nnew <- sum((x[ ,t-1] == 1)*(x[ ,t] == 2))
    N <- sum(x[ ,t-1] == 2)
    rho[t-1] <- Nnew / N
  }
  return(rho)
}))

par(mfrow = c(1,1))
plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "rho")
axis(1, c(1:nYears),labels = years)
myCol <- "firebrick3"
for(t in 1:(nYears-1)){
  plot.violins(list(posterior.rho[ ,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t

### ====       2.2.5.w ==== 
par(mfrow = c(1,1))
plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "w")
axis(1, c(1:nYears),labels = years)
for(t in 1:(nYears-1)){
  plot.violins(list(myResults$sims.list$w[ ,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t

params <- dimnames(nimOutput[[1]])[[2]][grep("w",dimnames(nimOutput[[1]])[[2]])]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}

### ====       2.2.6.phi ====  
par(mfrow = c(1,1))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "phi")
axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
myCol <- c("firebrick3")
for(t in 1:(nYears-1)){
  plot.violins(list(myResults$sims.list$phi[ ,t]),
               x = t,
               at = t,
               violin.width = 0.2,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t

params <- dimnames(nimOutput[[1]])[[2]][grep("phi",dimnames(nimOutput[[1]])[[2]])]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}

### ====       2.2.7.p0 ====  
## by country and trap-response
par(mfrow = c(1,2))
myDev <- c(-0.3, 0.3)
myCol <- c("blue4", "yellow3")

COUNTIES_AGGREGATEDSubsetsimp <- gSimplify(COUNTIES_AGGREGATEDSubset,tol = 500,topologyPreserve = T)
COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique

for(c in 1:dim(myResults$sims.list$p0)[2]){
  plot(myStudyArea)
  plot(COUNTIES_AGGREGATEDSubset[COUNTIES_AGGREGATEDSubsetsimp$idunique %in% c, ], add=T, col="red")
  
  
  plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
  axis(1, at = 1:(nYears), labels = years[1:(nYears)])
  
  for(t in 1:nYears){
    plot.violins(list(myResults$sims.list$p0[ ,c,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = "red",
                 add = T,
                 alpha = 0.3,
                 border.col = "red")
    # }#g
  }#t
}#c
# legend(x = 0, y = 0.25, legend = c("Norway", "Sweden"), fill= myCol)

params <- dimnames(nimOutput[[1]])[[2]][grep("p0", dimnames(nimOutput[[1]])[[2]])[-1]]
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}



### ====       2.2.8.the rest ====
params <- c("sigma", "dispSigma", "betaDens", "betaResponse", "betaCovs[1]", "betaCovs[2]", "betaCovs[3]")
for(i in 1:length(params)){
  PlotJagsParams(jags.samples = nimOutput, params = params[i])
}#i
dev.off()

### ==== 3. PLOT DENSITY MAPS ====
gc()
habbR <- raster::disaggregate(myHabitat$habitat.r, fact = 1)
habbR[habbR > 0] <- 1
plot(habbR)

## CALCULATE DENSITY RASTERS
Densi <- list()
for(t in 1:nYears){
  Densi[[t]] <- EstimateN_v3( habRaster = habbR, 
                              posterior.sxy = sxyScaled$data.scaled.xy[,,,t],
                              posterior.z = myResults_sxy$sims.list$z[,,t],
                              alive.states = 2,
                              return.all = FALSE,
                              regionEstimates = FALSE)
  print(t)
}#t

## PLOT
max.cuts <- max(unlist(lapply(Densi, function(x) max(x$PosteriorsCellsMean))))
cuts <- round(seq(0, max.cuts, length.out = 101) , digits = 2)
pal <- rev(terrain.colors(length(cuts)))
habdens <- habbR

pdf(file = file.path( myVars$WD, myVars$modelName,
                      paste(myVars$modelName, "_DENSITY_20km.pdf", sep = "")))
for(t in 1:nYears){
  habdens[] <- Densi[[t]]$PosteriorsCellsMean
  
  plot( habdens, breaks = cuts, col = pal,
        main = paste( years[t], "; N =", round(sum(Densi[[t]]$PosteriorsCellsMean[], na.rm = T))),
        legend = F, axes = F)   
  plot(myStudyArea, add = T, border = grey(0.5))
  plot( habdens, legend.only = TRUE, breaks = cuts, col = pal,
        legend.width = 2,
        axis.args = list(at = round(seq(0, max.cuts, length.out = 6), digits = 1),
                         labels = round(seq(0, max.cuts, length.out = 6), digits = 1),
                         cex.axis = 0.6),
        legend.args = list(text = 'Density', side=4, font=2, line=2.5, cex=0.8))
}
dev.off()

### ==== 4. PLOT DETECTIONS ====
pdf(file = file.path( myVars$WD, myVars$modelName,
                      paste(myVars$modelName, "_DETECTIONS.pdf", sep = "")))
myLayOut.mx <- cbind(c(1), c(1))
myLayOut <- layout(myLayOut.mx, width = c(1), heights = c(1))
#layout.show(myLayOut)

for(t in 1:nYears){
  ## PLOT STUDY AREA
  par(mar = c(1,1,1,1))
  plot(myBufferedArea, main = years[t], col = rgb(34/250, 139/250, 34/250, alpha = 0))
  plot(GLOBALMAP, col = "gray80", add = TRUE)
  plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
  plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
  
  ## PLOT DETECTIONS
  points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ], pch = 3, col = "darkred", cex = 0.5)
  points(myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ], pch = 3, col = "darkblue", cex = 0.5)
  
  
  ## ADD NUMBER OF INDIVIDUALS DETECTED
  graphics::text(x = 190000, y = 7880000, cex = 1,
                 labels = paste(length(unique(myFilteredData.sp$alive$Id[myFilteredData.sp$alive$Year == years[t]])), "Individuals"))
  
  
  ## ADD NUMBER OF DETECTIONS
  graphics::text(x = 190000, y = 7820000, cex = 1,
                 labels = paste(dim(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ])[1], "NGS samples"))
  
  
  ## ADD NUMBER OF DEAD RECOVERIES
  graphics::text(x = 190000, y = 7760000, cex = 1,
                 labels = paste(dim(myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ])[1], "Dead Recoveries"))
  
}#t
dev.off()

## -----------------------------------------------------------------------------------------------