#' @title State Matrix (z) Simulation function based on an Individual Based processes
#'
#' @description
#' \code{SimulateZIndividualBased} returns a \code{matrix} object with the individual states for each time step as well as its AC locations.
#' The states are simulated based on a post reproduction "matrix". This means that survey (the states) are determined just after the reproduction 
#' (i.eindividuals need to survive before reproducing )
#'
#' @param Initial.N A \code{vector} object containing the number of individuals to enter in each state
#' @param PHI A \code{character} matrix representing the state transtion (e.g. phi, r,...) 
#' @param REPRO A \code{character} matrix representing the transition of the reproduction 
#' @param FEC A \code{vector} representing the fecundity rate. A poisson process is applied to this rate to determine the number of offspring. 
#' @param n.occasions A \code{numeric} object to define the number occasions to be simulated.
#' @param Alive.states A \code{vector} with the different states considered alive. 
#' @param intercepts  A \code{list}, where for each of the matrices (PHI and REPRO), 
#' there is a list of  intercepts to be used for each parameter of the matrices (e.g. list(PHI=list( phi =0.5, r = 0.2), REPRO=list( pRepro=0.8)).
#' Parameter names should match the parameters entered in the matrices. Intercept are expressed on a probability scale.
#' @param id.betas A \code{list}, where for each of the matrices (PHI and REPRO),
#'  there is a list of betas to be used for each spatial covariate that where covariate names matches the \code{stack.cov.r} 
#'  (e.g.list( PHI = list( phi = list(longitude = +1, latitude  = +2 ), r = list(longitude = -1, latitude  = +0.5))). Betas are on the logit scale. 
#' @param stack.cov.r A \code{rasterstack} with the spatial covariates to be used to extract individual covariates for each individual and be used as a covariate.
#' @param id.betas A \code{list} of \code{vector} with beta coefficients on the log scale for each parameter. The interecept used is the one entered in the \code{intercepts}. The parameters should match the names entered in the state.trans object
#' @param habitat.r A \code{raster} with extend and resolution matching the stack.cov.r that represents available habitat (0/1) 
#' @param habitat.sp A \code{SpatialPoints} with the coordinates the raster cells.  
#' 
#' @param AC.method A \code{character} defining the method to be used by \code{SimulateMultiYearACs}. 
#' @param Hyper.sigma see \code{SimulateMultiYearACs} function 
#' @param covs.AC.dens see \code{SimulateMultiYearACs} function 
#' @param betas.AC.dens see \code{SimulateMultiYearACs} function 
#' @param habitat.poly if \code{AC.method="random"} provide a polygon from where the ACs will be  sampled.
#' @param plot.check A \code{logical} to display z image and covariates effect

#' @return z; A \code{matrix}  z object with individual state-histories z.
#'         AC.arr; A \code{array} with the AC coordinates of each ID at each occasion
#'         individuals.arr; A \code{individuals.arr} with the all information stored for each ID at each occasions 
#'         PHI.REALIZED; A \code{array} with the realised PHI matrix at each time step
#'         FECO.REALIZED; A \code{matrix} with the realised FECO vector at each time step
#'         REPRO.REALIZED; A \code{array} with the realised REPRO matrix at each time step   
#'         NNew; A \code{vector} with the number of new individuals at each time step         
#'         id.covs.phi; A \code{array} with covariates values used for each individual for the phi matrix              
#'         id.covs.repro; A \code{array} with covariates values used for each individual for the repro matrix             
#'         
#' @examples
#' 
#'  CREATE A RASTER HABITAT
#' habitat.r <- raster(matrix(1, nrow=10, ncol=10))
#' plot(habitat.r)
#' 
#' #==== 1. EXAMPLE WITHOUT  COVARIATES ====
#' 
#' # CREATE THE PHI MATRIX (survival and dead recovery)
#' # STATE 1 = ALIVE
#' #       2 = RECENTLY DEAD
#' #       3 = DEAD
#' PHI <-   matrix(c("phi"          ,    "0"   ,     "0",
#'                   "(1-phi)*r"    ,     "0"  ,     "0",
#'                   "(1-phi)*(1-r)",    "1"   ,    "1"), nrow=3, byrow=TRUE)
#' 
#' # CREATE THE REPRO MATRIX (FOR EACH STATES)
#' REPRO <- matrix(c("pRepro", "0"   ,    "0",
#'                   "0"     , "0"   ,    "0",
#'                   "0"     , "0"   ,    "0"), nrow=3,byrow=TRUE)
#' 
#' # CREATE THE FECUNDITY VECTOR (FOR EACH STATES)
#' FEC <- c(2, 0, 0)
#' 
#' # INITIAL POPULATION VECTOR  
#' Initial.N <- c(50, 0, 0)
#' 
#' ## SIMULATE IBM
#' sim.states <- SimulateZIndividualBased(  n.occasions = 4
#'                                          ,Initial.N = Initial.N
#'                                          ,habitat.sp = habitat.sp
#'                                          ,PHI=PHI
#'                                          ,REPRO=REPRO
#'                                          ,FEC=FEC
#'                                          ,AC.method="hyper"
#'                                          ,Alive.states=c(1)
#'                                          ,Hyper.sigma=0.01
#'                                          ,plot.check=TRUE
#'                                          ,covs.AC.dens=NULL
#'                                          ,betas.AC.dens=NULL
#'                                          ,habitat.r = habitat.r
#'                                          ,id.betas = list( PHI=NULL
#'                                                            , REPRO= NULL )
#'                                          ,stack.cov.r = list(  PHI=NULL
#'                                                                , REPRO=NULL)
#'                                          ,intercepts = list( PHI=list(  phi = 0.8
#'                                                                         , r = 0.2)
#'                                                              , REPRO=list( pRepro=0.5))
#' )
#' 
#' #==== 2. EXAMPLE WITH COVARIATES ====
#' # create covariates 
#' coord <- coordinates(habitat.r)
#' habitat.sp <- data.frame(coord)
#' coordinates(habitat.sp) <- habitat.sp
#' 
#' r.1 <- rasterFromXYZ(cbind(coord, coord[,1]))
#' r.2<- rasterFromXYZ(cbind(coord, coord[,2]))
#' r <- stack(r.1,r.2)
#' r <- scale(r)
#' names(r) <- c("longitude", "latitude")
#' plot(r)
#' 
#' 
#' 
#' ## SIMULATE IBM
#' 
#' sim.states <- SimulateZIndividualBased(  n.occasions = 4
#'                                          , Initial.N = Initial.N
#'                                          , habitat.sp = habitat.sp
#'                                          , PHI=PHI
#'                                          , REPRO=REPRO
#'                                          , FEC=FEC
#'                                          , Alive.states=c(1)
#'                                          , Hyper.sigma=0.01
#'                                          , plot.check=TRUE
#'                                          , covs.AC.dens=NULL
#'                                          , betas.AC.dens=NULL
#'                                          , habitat.r = habitat.r
#'                                          , id.betas = list( PHI=list( phi = list(  longitude = -1
#'                                                                                    , latitude  = 2 )
#'                                                                       , r = list(  longitude = 0.5
#'                                                                                    , latitude  = 1 ))
#'                                                             , REPRO= NULL )
#'                                          , stack.cov.r = list(  PHI=r
#'                                                                 , REPRO=NULL)
#'                                          , intercepts = list( PHI=list(  phi = 0.7
#'                                                                          , r = 0.2)
#'                                                               , REPRO=list( pRepro=0.5)))
#' 

SimulateZIndividualBased <- function(
                                      n.occasions = 10
                                      , 
                                      Initial.N = Initial.N
                                      , 
                                      habitat.sp = myHabitat.list$habitat.sp
                                      ,
                                      PHI=PHI
                                      ,
                                      REPRO=REPRO
                                      ,
                                      FEC=FEC
                                      ,
                                      Alive.states=c(1)
                                      ,
                                      AC.method="hyper"
                                      ,
                                      Hyper.sigma=0.2
                                      ,
                                      plot.check=TRUE
                                      ,
                                      covs.AC.dens=NULL
                                      ,
                                      betas.AC.dens=NULL
                                      ,
                                      habitat.r=myHabitat$habitat.r#IDCells.mx= myHabitat.list$
                                      ,
                                      id.betas= list(PHI=list( phi=list(layer.1=c(-2))
                                                               , r =list(layer.1= c(+2)))
                                                     ,REPRO=list( pRepro=list(layer.2=c(-2))))
                                      ,
                                      stack.cov.r = list( PHI=r.phi
                                                          ,REPRO=r.phi)
                                      ,
                                      intercepts= list(PHI=list( phi=0.8
                                                                 , r= 0.2)
                                                       ,REPRO=list( pRepro=0.8))
                                      ,
                                      habitat.poly=NULL
                                      ,
                                      fix.f=FALSE
                                      
){
  #require(popbio)
   
  ## -------------------------------------------------------------------------------------------
  ## ------ STEP 1 :  #SET UP THE THE NECESSARY OBJECTS AND INITIALIZE POP AT TIME STEP 1 ------ 
  ## -------------------------------------------------------------------------------------------
  # INITIALIZE POPULATION SIZE AT TIME STEP 1
  # KEEP TRACK OF # INDIVIDUALS 
  Nt <- NNew <- 0
  Nt[1] <- sum(Initial.N)
  NNew[1] <- sum(Initial.N)
  
  # CREATE A TABLE TO KEEP TRACK OF THE INDIVIDUALS 
  M <- list()
  M[[1]] <- matrix(NA, nrow=Nt[1], ncol=9)
  colnames(M[[1]]) <- c("ID","age", "z", "r","f", "Parent.ID","sx","sy","IDCell")
  M[[1]][,"ID"] <- 1 : Nt[1]
  M[[1]][,"age"] <- 1
  
  Z <- list() # THE STATE OF IDS 
  R <- list() # WHETHER THEY REPRODUCE OR NOT 
  B <- list() # NUMBER OF OFFSPRINGS 
  
  # UPDATE STATE STATUS
  M[[1]][,"z"] <- rep(c(1:length(Initial.N)), Initial.N )  
  Z[[1]] <- rep(c(1:length(Initial.N)), Initial.N ) 
  
  
  # ASSIGN AC LOCATIONS TO ALL IDS FROM FIRST OCCASION TO LAST OCCASIONS 
  if(AC.method=="random"){
     ACS <- list()
     ACS$AC.sp.list <- list()
     
     AC <- spsample( habitat.poly, Nt[1], type="random")
     for(tt in 1 :n.occasions){
        ACS$AC.sp.list[[tt]] <- AC 
     }
     
     }else{
     ACS <-  SimulateMultiYearACs ( N = Nt[1]
                                    , n.occasions = n.occasions
                                    , habitat.sp=habitat.sp
                                    , covs = covs.AC.dens
                                    , betas = betas.AC.dens
                                    , type = "Hyper.AC"
                                    , sigma=Hyper.sigma
                                    , initial.AC.sp = NULL
                                    , plot.check = FALSE)
  }
  # STORE ACS IN AN ARRAY 
  AC.arr <- array(NA, c(Nt, 2,n.occasions))
  dim(AC.arr)
  
  dimnames(AC.arr)[2] <- list(c("x","y"))
  for(t in 1:length(ACS$AC.sp.list)){
    AC.arr[,,t] <- coordinates(ACS$AC.sp.list[[t]] )
  }
  
  # STORE  SX SY COORDINATES IN THE TABLE 
  M[[1]][,c("sx", "sy")] <- coordinates(ACS$AC.sp.list[[t]])
  
  # STORE THE INDIVIDUAL COVARIATES 
  id.covs.repro <- list()
  id.covs.phi <- list()
  
  ind.fixed.effects.PHI <- list()
  ind.fixed.effects.REPRO <- list()
  ## ----------------------------------------
  ## ------ STEP 2 : FOR EACH OCCASION ------ 
  ## ----------------------------------------
  
  for(t in 2:n.occasions){
    # OBTAIN CELL ID IN WHICH AC IS LOCATED 
    for(i in 1:nrow(M[[t-1]])){
      M[[t-1]][i,c("IDCell")] <-  cellFromXY(habitat.r, M[[t-1]][i, c("sx","sy")] ) #IDCells.mx[trunc(M[[t-1]][i,c("sy")] )+1, trunc(M[[t-1]][i,c("sx")] )+1]
    }
    
     
     ## ====================================================
     ## ==== STEP 2.1 : INTEGRATE INDIVIDUAL COVARIATES ==== 
     ## ====================================================
       # ==== 2.1.1. PHI MATRIX ==== 
    # extract covariate values 
    if(!is.null(id.betas$PHI)){
      id.covs.phi[[t-1]] <-  data.frame(stack.cov.r$PHI[M[[t-1]][,c("IDCell")]])
      colnames(id.covs.phi[[t-1]]) <- names(stack.cov.r$PHI)
    }
     
     
    # A PHI MATRIX FOR EACH INDIVIDUAL
    OMEGA.PHI <- array(NA, c(dim(PHI)[1], dim(PHI)[1], nrow(M[[t-1]])))
    OMEGA.PHI[ ] <- PHI
    
    for(p in 1:length(intercepts$PHI)){
        if(!is.null(id.betas$PHI)){# if covariates 
      temp <- formula(paste("~", paste(names(id.betas$PHI[[p]]), collapse="+"), "-1", sep=""))
      df.id <- data.frame(id.covs.phi[[t-1]][,names(id.betas$PHI[[p]])])
      colnames(df.id) <- names(id.betas$PHI[[p]])
      
      Xmat <- model.matrix(temp, df.id)
      ind.fixed.effects <- inv.logit(Xmat %*% unlist(id.betas$PHI[[p]][names(id.betas$PHI[[p]])]) + logit(intercepts$PHI[[p]]))
      ind.fixed.effects.PHI[[names(intercepts$PHI)[p]]][[t-1]] <- ind.fixed.effects
        }else{#intercept only
        ind.fixed.effects <-  matrix(intercepts$PHI[[p]], nrow=nrow(M[[t-1]]))
        }
       
      # fill in the individual i matrix with the values   
      for(i in 1:nrow(M[[t-1]])){
        OMEGA.PHI[,,i] <- gsub(names(intercepts$PHI)[p]
                               , ind.fixed.effects[i,]
                               , OMEGA.PHI[,,i]) 
      }#i
    }#p
    
    OMEGA.PHI <- apply(OMEGA.PHI, c(1:3), function(x) eval(parse(text=x))  )
    
 
      # ==== 2.1.2. REPRO MATRIX  ====
    # extract covariate values 
    if(!is.null(id.betas$REPRO)){
      id.covs.repro[[t-1]] <-  data.frame(stack.cov.r$REPRO[M[[t-1]][,c("IDCell")]])
      colnames(id.covs.repro[[t-1]]) <- names(stack.cov.r$REPRO)
      
    }
    
    # A REPRO MATRIX FOR EACH INDIVIDUAL
    
    OMEGA.REPRO <- array(NA, c(dim(REPRO)[1], dim(REPRO)[1], nrow(M[[t-1]])))
    OMEGA.REPRO[ ] <- REPRO
    
    for(p in 1:length(intercepts$REPRO)){
      if(!is.null(id.betas$REPRO)){# if covariates 
      temp <- formula(paste("~", paste(names(id.betas$REPRO[[p]]), collapse="+"), "-1", sep=""))
      df.id <- data.frame(id.covs.repro[[t-1]][, names(id.betas$REPRO[[p]])])
      colnames(df.id) <- names(id.betas$REPRO[[p]])
      
      Xmat <- model.matrix(temp, df.id)
      ind.fixed.effects <- inv.logit(Xmat %*% unlist(id.betas$REPRO[[p]][names(id.betas$REPRO[[p]])]) + logit(intercepts$REPRO[[p]]))
      ind.fixed.effects.REPRO[[names(intercepts$REPRO)[p]]][[t-1]] <- ind.fixed.effects
      }else{#intercept only
      ind.fixed.effects <-  matrix(intercepts$REPRO[[p]], nrow=nrow(M[[t-1]]))
      }
      
      # fill in the individual i matrix with the values   
      for(i in 1:nrow(M[[t-1]])){
        OMEGA.REPRO[,,i] <- gsub(names(intercepts$REPRO)[p]
                                 , ind.fixed.effects[i,]
                                 , OMEGA.REPRO[,,i] ) 
        
      }#i
    }#p
    
    # plot(r.phi[[2]])
    # i=50
    # points(M[[t-1]][i,c("sy")] ~ M[[t-1]][i,c("sx")] )
    # OMEGA.REPRO[,,i]
    
    OMEGA.REPRO <- apply(OMEGA.REPRO, c(1:3), function(x) eval(parse(text=x))  )
    ## ====================================================
    ## ==== STEP 2.2 : POPULATION DYNAMICS PROCESS     ==== 
    ## ====================================================
    # LIST FOR INDIVIDUAL STATE Z 
    Z[[t]] <- rep(NA,length(Z[[t-1]]))
    # LIST FOR INDIVIDUAL REPRODUCTION STATE (0/1)
    R[[t]] <- rep(NA,length(Z[[t-1]]))
    # LIST FOR INDIVIDUAL NUMBER OF OFFSPRINGS (Poisson)
    B[[t]] <- rep(NA,length(Z[[t-1]]))
    
    # FOR EACH ID
    for(i in 1:length(Z[[t-1]])){
      # SURVIVAL PART 
      Z[[t]][i] <- which(rmultinom(1, 1, OMEGA.PHI[,Z[[t-1]][i],i])==1)
      
      # IF SURVIVED, REPRO?
      if(sum(Z[[t]][i] %in% Alive.states)==1){
        R[[t]][i] <- rbinom(1, 1, OMEGA.REPRO[1,Z[[t]][i],i])
      }else{
        R[[t]][i] <- 0
      }
      # IF REPRO, HOW MANY PUPS? 
      if(R[[t]][i]==1 ){
         if(fix.f==T){
            B[[t]][i] <- 1 
         }else{
        B[[t]][i] <- rpois(1, FEC[Z[[t]][i]])
         }
      }else{
        B[[t]][i] <- 0
      }
    }#i
    
    ## ======================================================
    ## ==== STEP 2.3 : UPDATE THE OBJECTS WITH NEW BORNS ==== 
    ## ======================================================
    # COUNT THE NEWBORN AND ALIVE IND
    NNew[t] <- sum(B[[t]])
    Nt[t] <-  NNew[t]+ sum(Z[[t]] %in% Alive.states)
    
    # ==== 2.3.1. GET THE OLD GUYS AND UPDATE THEIR STATES 
    m <- M[[t-1]]   
    # UPDATE THEIR SATE, R AND F STATES 
    m[, c("z","r","f")] <- cbind(Z[[t]],R[[t]],B[[t]])
    # UPDATE THEIR AGE
    m[,"age"] <- m[,"age"] + 1 
    
    # UPDATE THEIR ACS
      # COORDINATES
    m[,c("sx","sy")] <- AC.arr[,,t]
      # CELL ID
    for(ii in 1:nrow(m)){
      m[ii,c("IDCell")] <-  cellFromXY(habitat.r, m[ii,c("sx","sy")] ) #IDCells.mx[trunc(m[ii,c("sy")] )+1, trunc(m[ii,c("sx")] )+1]
    }#ii
    
    # ==== 2.3.2. GET THE NEW BORNS AND UPDATE THEIR STATES 
    m.new <- matrix(NA, nrow=NNew[t], ncol=9)
    colnames(m.new) <- c("ID","age", "z", "r","f", "Parent.ID","sx","sy","IDCell")
    m.new[,"ID"] <- (nrow(M[[t-1]])+1) : (NNew[t]+(nrow(M[[t-1]]))) 
    m.new[,"z"] <- 1
    m.new[,"age"] <- 1
    Z[[t]] <- c(Z[[t]], m.new[,"z"] )
    
    m.new[ ,"Parent.ID"] <- rep(m[,"ID"], B[[t]])
    
    # ASSIGN AC OF THE PARENTS TO THE NEW BORNS 
    m.new[,"sx"] <- m[m.new[,"Parent.ID"],"sx"]
    m.new[,"sy"] <- m[m.new[,"Parent.ID"],"sy"]
    
   
    
    ## ======================================================
    ## ==== STEP 2.4 : UPDATE THE OBJECTS WITH NEW BORNS ==== 
    ## ======================================================
    # SIMULATE MULTI-YEAR AC FOR THE NEW BORNS 
    # THEY GET THEIR PARENTS AC AS AN HYPER AC.
    
    if(AC.method=="random"){
       ACS.1 <- list()
       ACS.1$AC.sp.list <- list()
       
       AC <- spsample( habitat.poly, NNew[t], type="random")
       for(tt in 1 :n.occasions){
          ACS.1$AC.sp.list[[tt]] <- AC 
       }
       
       #update the first year AC location. so it is define randomly
       sxy.newborn.sp <- data.frame(m.new[,c("sx","sy")])
       sxy.newborn.sp[,1:2] <- coordinates(ACS.1$AC.sp.list[[1]])
       m.new[,c("sx","sy")] <- coordinates(ACS.1$AC.sp.list[[1]])
       coordinates(sxy.newborn.sp) <- sxy.newborn.sp
       proj4string(sxy.newborn.sp) <- CRS(proj4string(habitat.sp))
       
       AC.arr1 <- array(NA, c(NNew[t], 2, n.occasions))
       AC.arr1[,,t] <- m.new[,c("sx","sy")]
       
       
       if(t!=n.occasions){
          for(t1 in (t+1):n.occasions){
             AC.arr1[,,t1] <- coordinates(ACS.1$AC.sp.list[[t1-(t)]] )
          }#t1
       }
       

       }else{
    
      sxy.newborn.sp <- data.frame(m.new[,c("sx","sy")])
      coordinates(sxy.newborn.sp) <- sxy.newborn.sp
      proj4string(sxy.newborn.sp) <- CRS(proj4string(habitat.sp))
      AC.arr1 <- array(NA, c(NNew[t], 2, n.occasions))
      AC.arr1[,,t] <- m.new[,c("sx","sy")]
      
      if(t!=n.occasions){## dont simulate ac for the last occasion 
         
      ACS.1 <-  SimulateMultiYearACs ( N = NNew[t]
                                       , n.occasions = n.occasions- (t)
                                       , habitat.sp=habitat.sp
                                       , covs = covs.AC.dens
                                       , betas = betas.AC.dens
                                       , type = "Hyper.AC"
                                       , sigma = Hyper.sigma
                                       , initial.AC.sp = sxy.newborn.sp
                                       , plot.check = FALSE)
      
         # FILL AC COORDINATES IN AN ARRAY
         for(t1 in (t+1):n.occasions){
            AC.arr1[,,t1] <- coordinates(ACS.1$AC.sp.list[[t1-(t)]] )
         }#t1
      
      }
      
    }
   
       
    # UPDATE CELL ID 
    for(ii in 1:nrow(m.new)){
       m.new[ii,c("IDCell")] <- cellFromXY(habitat.r, m.new[ii,c("sx","sy")] ) # IDCells.mx[trunc(m.new[ii,c("sy")] )+1, trunc(m.new[ii,c("sx")] )+1]
    }#ii
    
    # MERGE NEW BORNS AND OLD 
    M[[t]] <- rbind(m, m.new)
    
     # BIND AC OLD AND NEW BORNS 
      AC.arr <-  abind(AC.arr, AC.arr1, along = 1)
  }#t
  
  
  #GET THE FINAL TABLE IN A ARRAY
  arr <- array(NA, c(dim(M[[n.occasions]])[1], ncol(M[[1]]), n.occasions ))
  for(t in 1:length(M)){
    arr[1:dim(M[[t]])[1],,t]<- M[[t]]
  }#t
  
  dimnames(arr)[2] <- list(c("ID","age", "z", "r","f", "Parent.ID","sx","sy","IDCell"))
  
  
  # EXTRACT THE STATE OF INDIVIDUALS AT EVERY TIME STEP 
  z <- arr[,"z",]
  z[is.na(z)] <- 0
  
  # EXTRACT THE REPRO OF INDIVIDUALS AT EVERY TIME STEP 
  repro <- arr[ ,"r", ]
  feco <- arr[ ,"f", ]
  
  
  
 # OBTAIN REALISED TRANSITION
  n.states <- 1:dim(PHI)[1]
   
   # FOR PHI
  PHI.REALIZED <- array(NA, c(dim(PHI),(dim(z)[2]-1) ))
  
  for(t in 2:dim(z)[2]){
   for(s in 1: dim(PHI)[1]){
     n.states1 <- n.states[!( PHI[,s] %in% c("0"))]
     for( s1 in 1:length(n.states1)){
        which. <-  which(z[,t-1]==n.states[s])
        which.t <-  which(z[which.,t]==n.states1[s1])
        PHI.REALIZED[n.states1[s1],n.states[s], t-1] <- length(which.t)/length(which.)
      }
    }
  }
  
  # FOR feco

  FECO.REALIZED <- array(NA, c(length(FEC),(dim(feco)[2]) ))
  
  for(t in 1:dim(feco)[2]){
        n.states1 <- n.states[!( FEC %in% c("0"))]
        for( s1 in 1:length(n.states1)){
           which.repro <-  which(repro[,t]==1)# you have to be alive to reproduce 
           FECO.REALIZED[s1,t] <-  sum(feco[which.repro,t])/length(which.repro)
        }  
     }

  
 
  # FOR REPRO
 
  REPRO.REALIZED <- array(NA, c(dim(REPRO),(dim(repro)[2]) ))
  
  for(t in 2:dim(repro)[2]){
     for(s in 1: dim(REPRO)[1]){
        n.states1 <- n.states[!( REPRO[,s] %in% c("0"))]
        for( s1 in 1:length(n.states1)){
           which.alive <-  which(z[,t]==n.states[s])
           REPRO.REALIZED[n.states1[s1], n.states[s], t] <-  mean(repro[which.alive, t],na.rm=T )
        }  
     }
  }
  
  
   
  ## ==============================
  ## ==== STEP 3 : CHECK PLOTS ==== 
  ## ==============================
  if(plot.check){
    par(mfrow=c(1,2), mar=c(3,3,3,3))
    # PLOT POP TRAJECTORY  
    GetPopCompo(z, cex = 0.7, legend.placement = "topleft")
    
   
    ## COMPARE IT TO MATRIX PROJECTION
    # mat <- apply(OMEGA.PHI,c(1,2), mean)
    # mat[1,] <- mat[1,] * FEC + apply(OMEGA.REPRO,c(1,2), mean)[1,]
    # p <- pop.projection(mat, Initial.N, n.occasions)
    # points(p$stage.vectors[Alive.states,], type="l", lty=2, col="red")

    
    #PLOT ALL YEARLY ACS FROM ALIVE IND
    par( mar=c(0,0,0,0))
    col <- adjustcolor(rainbow(n.occasions), alpha.f = 0.4)
    plot(coordinates(habitat.sp)[,2]~coordinates(habitat.sp)[,1], col=grey(0.9))
    legend("topright", legend =1:n.occasions, col=col, pch=16, title = "occasions", bty = "n"  )
    
    for(t in 1:n.occasions){
      m.plot <- M[[t]] 
      alive <- m.plot[m.plot[,"z"]==1,]
      points(alive[,c("sy")]~alive[,c("sx")], col=col[t], pch=16)
    }
    
    # PLOT REPRO AND PHI COVARIATE EFFECTS 
    #id.betas <- id.betas[[1]]
    for(what in 1:length(id.betas)){
    
    if(!is.null(id.betas[[what]])){
       
       for(p in 1:length(id.betas[[what]])){
          par(mfrow=c(1, length(id.betas[[what]][[p]])+1), mar=c(5,5,3,3))
          
          for(ncov in 1: length(id.betas[[what]][[p]]) ){
             cov.values <- do.call(rbind, id.covs.phi)
             myX <- seq(min(cov.values[,names(id.betas$PHI[[p]])[ncov]]), max(cov.values[,names(id.betas[[what]][[p]])[ncov]]), length.out = 1000)
             myY <- inv.logit(logit(intercepts[[what]][[p]]) + id.betas[[what]][[p]][[ncov]] * myX)
             
             color <- rev(terrain.colors(length(myX)))
                 
             plot(myX, myY, xlab = names(id.betas[[what]][[p]])[ncov], ylab = names(id.betas[[what]])[p],
                  ylim = c(0,max(myY)), col =color )
          }#ncov
       
       par(mar=c(2,0,2,2))
       #plot(habitat.sp, col=grey(0.9))
       
     
       r.phi <- stack.cov.r[[what]] 
       r.phi1 <- r.phi[[1]]
       
       ind.fixed.effects.R <- inv.logit(matrix(r.phi[[names(id.betas[[what]][[p]])]], nrow =ncell(r.phi) ) %*% unlist(id.betas[[what]][[p]][names(id.betas[[what]][[p]])]) + logit(intercepts[[what]][[p]]))
       r.phi1[] <- ind.fixed.effects.R[,1]
       
       col.breaks <- seq(min(ind.fixed.effects.R[,1],na.rm = T), max(ind.fixed.effects.R[,1],na.rm = T), length.out = 1000)
       colo <- rev(terrain.colors(length(col.breaks)))
       
       plot(r.phi1, main=paste("TRUE", names(id.betas[[what]])[p]), breaks=col.breaks, col=colo
            , axis.args=list(at=round(seq(min(ind.fixed.effects.R, na.rm = T), max(ind.fixed.effects.R, na.rm = T), length.out = 7), digits = 3),
                           labels=round(seq(min(ind.fixed.effects.R, na.rm = T), max(ind.fixed.effects.R, na.rm = T), length.out = 7), digits = 3), 
                           cex.axis=0.6))
       
       for(t in 1:(n.occasions-1)){
          m.plot <- M[[t]] 
          ID <- which(m.plot[,"z"]==1)
          if(names(id.betas)[what]=="REPRO") {
             id.ind.fixed.effects <-  ind.fixed.effects.REPRO[[names(id.betas[[what]])[p]]][[t]][ID,]
               }else{
                  id.ind.fixed.effects <-  ind.fixed.effects.PHI[[names(id.betas[[what]])[p]]][[t]][ID,]
                  
               }
     
          points(m.plot[ID,c("sy")] ~ m.plot[ID,c("sx")], pch=16, col="black", cex=0.8)
          
          if(sum(id.ind.fixed.effects != r.phi1[m.plot[ID, c("IDCell")]]) != 0){ print("WARNINGS! TRUE SURVIVAL AND TRUE SURVIVAL ASSIGNED TO IDS ARE NOT EQUALS")}
          }#t
      }#p
    }
}#what
  
    # PLOT AVERAGE REALISED TRANSITION 
    par(mfrow=c(1,3), mar=c(5,5,5,2))
      #PHI
    plot(1, type="n", xlab="Occasions", ylab="", main="Realised values; PHI matrix", xlim=c(0, dim(PHI.REALIZED)[3]), ylim=c(0, 1))
    
    trans.phi <- which( PHI != c("0"))
    trans.phi.coord <- which( PHI != c("0"),arr.ind = T)
    
    m <- matrix(NA,nrow=length(trans.phi), ncol=dim(PHI.REALIZED)[3])
    for( t in 1:dim(PHI.REALIZED)[3]){
       m[,t]<-  PHI.REALIZED[,,t][trans.phi]
    }
    col <- adjustcolor(rainbow(length(trans.phi)),alpha.f = 0.5)
    for( i in 1:length(trans.phi)){
       x <- c(1:dim(PHI.REALIZED)[3]) + jitter(0, factor =8)
         points(m[i,] ~ x, pch=16, col= col[i])
         points(m[i,] ~ x, pch=21, bg= col[i], type="b")
    }
    
    legend("topleft", col=col, legend = paste("from", trans.phi.coord[,2], "to", trans.phi.coord[,1]), pch=rep(16, length(trans.phi)), cex=0.7)
    
    
    #REPRO
    plot(1, type="n", xlab="Occasions", ylab="", main="Realised values; REPRO matrix", xlim=c(0, dim(REPRO.REALIZED)[3]), ylim=c(0, 1))
    
    trans.repro <- which( REPRO != c("0"))
    trans.repro.coord <- which( REPRO != c("0"),arr.ind = T)
    
    m <- matrix(NA,nrow=length(trans.repro), ncol=dim(REPRO.REALIZED)[3])
    for( t in 1:dim(REPRO.REALIZED)[3]){
       m[,t]<-  REPRO.REALIZED[,,t][trans.repro]
    }
    col <- adjustcolor(rainbow(length(trans.repro)),alpha.f = 0.5)
    for( i in 1:length(trans.repro)){
       x <- c(1:dim(REPRO.REALIZED)[3]) + jitter(0, factor =8)
       points(m[i,] ~ x, pch=16, col= col[i])
       points(m[i,]~ x, pch=21, bg= col[i], type="b")
    }
    
    legend("topleft", col=col, legend = paste("from", trans.repro.coord[,2], "to", trans.repro.coord[,1]), pch=rep(16, length(trans.repro)), cex=0.7)
    
    #F
    plot(1, type="n", xlab="Occasions", ylab="", main="Realised values; FEC (Per capita/fec) ", 
         xlim=c(0, dim(FECO.REALIZED)[2]), ylim=c(0, max(FECO.REALIZED, na.rm=T)))
    
    trans.fec <- which( FEC != c("0"))
    

    col <- adjustcolor(rainbow(length(trans.fec)),alpha.f = 0.5)
    for( i in 1:length(trans.fec)){
       x <- c(1:dim(FECO.REALIZED)[2]) + jitter(0, factor =8)
       points(FECO.REALIZED[i,] ~ x, pch=16, col= col[i])
       points(FECO.REALIZED[i,]~ x, pch=21, bg= col[i], type="b")
    }
    
    legend("topleft", col=col, legend = paste("from", trans.fec), pch=rep(16, length(trans.fec)), cex=0.7)
    
    
    
    
    
  }#plot.check
    
  individuals.arr <- arr
  return(list(  z = z
              , AC.arr = AC.arr
              , individuals.arr = arr
              , PHI.REALIZED = PHI.REALIZED
              , FECO.REALIZED = FECO.REALIZED
              , REPRO.REALIZED = REPRO.REALIZED
              , NNew = NNew
              , id.covs.phi=id.covs.phi
              , id.covs.repro=id.covs.repro))
  
}
  