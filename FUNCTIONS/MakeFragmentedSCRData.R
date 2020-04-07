#' @title Function to generate specific elements (inputs) needed to run a "Restricted" SCR model.

#'
#' @description
#' \code{MakeFragmentedSCRData} returns a list object with \code{} 
#' @param y \code{Numeric } detection history matrix, with individuals in rows 
#' (first detected individuals, then augmented individuals) and detectors in columns.
#' @param detector.xy \code{Numeric} Matrix with scaled detector coordinates (x and y in columns)
#' @param habitat.xy \code{Numeric} Matrix with habitat raster cell centroid coordinates (x and y in columns)
#' @param habitat.mx \code{Matrix} with of the habitat with 1: habitat and 0 : non-habitat
#' @param window.AC \code{Numeric} size of the sides of the square extent window that defines the boundaries of the individual AC regions
#' @param window.detectors \code{Numeric} size of the sides of the square extent window that defines the boundaries of the individual detectors regions
#' @param n.layers \code{Numeric} with the number of layers of augmented individuals to be added. 
#' @param plot.check A \code{logical} for whether (\code{TRUE}) or not (\code{FALSE}) plots are to be generated during simulations.
#' @param observable.states \code{Numeric} with the states that are observable in the y matrix (i.e. alive)
#' @param centroid.method \code{string} whether individual centroids of detected individuals should be defined annually ("yearly") or kept constant ("overall")
#' @param offset A \code{logical} if layers the exact same layer of augmented individuals should be added on top of each other (FALSE), or whether they should be added with an "offset" so they dont show perfect overlap 
#' @examples

MakeFragmentedSCRData  <- function( y
                                    , detector.xy
                                    , habitat.xy
                                    , habitat.mx
                                    , IDCells.mx
                                    , searched.mx 
                                    , window.AC = 2
                                    , window.detectors= 4
                                    , observable.states = 1
                                    , n.layers = 30
                                    , centroid.method = "overall"
                                    , plot.check = TRUE
                                    , augment.searched = FALSE
                                    , min.augment.layer=2
                                    , offset=FALSE)
{
  # ---- STEP 1: CHECK THE DIMENSIONS OF THE DATA & DECLARE VARIABLES ----- 
  if((window.AC> window.detectors)==1){  print("WARNINGS!!! Window.AC is larger than the window of the detectors")}
   
  if((window.AC %% 2)==1){ 
      window.AC <- window.AC+1 
   print(paste("WARNINGS!!! Window.AC is not an even value and has been rounded to", window.AC, " so it matches the habitat"))
   }
   
  if(dim(y)[2] != dim(detector.xy)[1]){stop(" The number of detectors does not match the detection array dimensions ! ")}

  if(augment.searched==TRUE){
     n.layers1 <- n.layers
     n.layers <- 1}
   
  n.years <- ifelse(length(dim(y))>=3, dim(y)[3], 1)
  n.individuals.detected <- dim(y)[1]
  n.detectors <- dim(detector.xy)[1]
  n.habitat <- dim(habitat.xy)[1]
  
  if(length(dim(y)) == 2){y <- array(y, c(n.individuals.detected, n.detectors, n.years))}
  if(length(dim(detector.xy)) == 2){ detector.xy.arr <- array(NA, c(n.detectors, 2, n.years))
                                       for(t in 1:n.years){
                                          detector.xy.arr[,,t] <- detector.xy
                                       }    
                                     detector.xy <- detector.xy.arr
                                    }   
  
  if(length(dim(habitat.xy)) == 2){ habitat.xy.arr <- array(habitat.xy, c(n.habitat, 2, n.years))
                                      for(t in 1:n.years){
                                         habitat.xy.arr[,,t] <- habitat.xy
                                      }    
                                     habitat.xy <- habitat.xy.arr
                                    }
  
  # ---- STEP 2: IDENTIFY INDIVIDUAL CENTROIDS FOR DETECTED INDIVIDUALS -----
  id.detectors <- list()
  xy.center.detected <- array(NA, c(n.individuals.detected, 2, n.years))
  
  # Method 1 : Overall Centroid
  if(centroid.method == "overall"){
    for(t in 1:n.years){	
      id.detectors[[t]] <- lapply(1:n.individuals.detected, function(x){which(y[x, , ] %in% observable.states)})
      xy.center.detected[ , ,t] <- do.call(rbind, lapply(id.detectors[[t]], function(x){
        x <- rbind(detector.xy[x, ,t], c(NA,NA))
        apply(x, 2, function(y) mean(y, na.rm=TRUE))}))  
    }#t
  }#if
  
  # Method 2 : Yearly Centroids
  if(centroid.method == "yearly"){
    for(t in 1:n.years){	
      id.detectors[[t]] <- lapply(1:n.individuals.detected, function(x){
        if(any(y[x, ,t] %in% observable.states)){
          which(y[x, ,t] %in% observable.states)
        }
        else{
          which(y[x, , ] %in% observable.states)
        }})
      xy.center.detected[ , ,t] <- do.call(rbind,lapply(id.detectors[[t]], function(x){
        x <- rbind(detector.xy[x, ,t], c(NA,NA))
        apply(x, 2, function(y)mean(y, na.rm=TRUE))}))         
    }#t
  }#if
  
  if(!centroid.method %in% c("yearly","overall"))stop(" ERROR : The method to calculate individual centroids must be one of /overall/ or /yearly/. ")
  
  
  
  # ---- STEP 3: IDENTIFY INDIVIDUAL CENTROIDS FOR AUGMENTED INDIVIDUALS -----
  # no offset
  if(offset!=TRUE){
  x.seq <- seq(round(min(habitat.xy[,1,])) - window.AC*2, round(max(habitat.xy[,1,])) + window.AC*2, by = window.AC)
  y.seq <- seq(round(min(habitat.xy[,2,])) - window.AC*2, round(max(habitat.xy[,2,])) + window.AC*2, by = window.AC)
  aug0.xy <- as.data.frame(expand.grid(x.seq, y.seq))
  
  aug.xy <- do.call(rbind, lapply(1:n.layers, function(x) aug0.xy ))
  
   }else{
  #offset
  if(window.AC < n.layers){## IF the window size is larger than the number of layers, add layers of layers.
    vec <- c(1:n.layers)
    split.vec <- split(vec, ceiling(vec/window.AC))
    
    aug.xy <- list()
    for( la in 1:length(split.vec)){
      shift.size <- round(window.AC/length(split.vec[[la]]))# round to fit the habitat 
      x.seq <- seq(round(min(habitat.xy[,1,])) - window.AC*2, round(max(habitat.xy[,1,])) + window.AC*2, by = window.AC)
      y.seq <- seq(round(min(habitat.xy[,2,])) - window.AC*2, round(max(habitat.xy[,2,])) + window.AC*2, by = window.AC)
      aug0.xy <- as.data.frame(expand.grid(x.seq, y.seq))
      aug.xy[[la]] <- do.call(rbind, lapply(1:n.layers, function(x) aug0.xy + (x-1)*shift.size))
    }
    
    aug.xy <- do.call(rbind, aug.xy)
    
  }else{
    shift.size <- round(window.AC/n.layers)# round to fit the habitat 
    x.seq <- seq(round(min(habitat.xy[,1,])) - window.AC*2, round(max(habitat.xy[,1,])) + window.AC*2, by = window.AC)
    y.seq <- seq(round(min(habitat.xy[,2,])) - window.AC*2, round(max(habitat.xy[,2,])) + window.AC*2, by = window.AC)
    aug0.xy <- as.data.frame(expand.grid(x.seq, y.seq))
    aug.xy <- do.call(rbind, lapply(1:n.layers, function(x) aug0.xy + (x-1)*shift.size))
     }
   }  
  
  n.individuals.augmented <- dim(aug.xy)[1]
  xy.center.augmented <- array(NA, c(n.individuals.augmented, 2, n.years))
  xy.center.augmented[ ,1, ] <- aug.xy[ ,1]
  xy.center.augmented[ ,2, ] <- aug.xy[ ,2]
  
  # ---- STEP 4: DELINEATE MOVING WINDOW FOR AC (MIN/MAX OF STATE SPACE X AND Y) ----- 
  n.individuals <- n.individuals.detected + n.individuals.augmented
  xy.center <- abind(round(xy.center.detected), xy.center.augmented, along = 1)
  MovingWindows.AC.xy  <- abind(xy.center, xy.center - window.AC/2, xy.center + window.AC/2, along = 2)
  
  dimnames(MovingWindows.AC.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  
  xy.bounds.AC <- array(NA, c(n.individuals, 2, 2, n.years))
  xy.bounds.AC[ , ,1, ] <- MovingWindows.AC.xy[,c("lower.x", "lower.y"),]
  xy.bounds.AC[ , ,2, ] <- MovingWindows.AC.xy[,c("upper.x", "upper.y"),]
  
  # ---- STEP 5: DELINEATE MOVING WINDOW FOR DETECTORS  ----- 
  MovingWindows.detectors.xy <- MovingWindows.AC.xy
  diff.windows <- (window.detectors-window.AC)/2
  for(t in 1:n.years){
     # extend moving windows for detectors 
     MovingWindows.detectors.xy[,c("lower.x","lower.y"),t] <- MovingWindows.detectors.xy[,c("lower.x","lower.y"),t] -  diff.windows    
     MovingWindows.detectors.xy[,c("upper.x","upper.y"),t] <- MovingWindows.detectors.xy[,c("upper.x","upper.y"),t] +  diff.windows 
  }
  
  for(t in 1:n.years){
     # cut to habitat extent
     xmin <- min(floor(habitat.xy[,1,t]))
     ymin <- min(floor(habitat.xy[,2,t]))
     xmax <- max(ceiling(habitat.xy[,1,t]))
     ymax <- max(ceiling(habitat.xy[,2,t]))
     
     MovingWindows.detectors.xy[MovingWindows.detectors.xy[ ,"lower.x",t] < xmin, "lower.x", t] <- xmin
     MovingWindows.detectors.xy[MovingWindows.detectors.xy[ ,"lower.y",t] < ymin, "lower.y", t] <- ymin
     MovingWindows.detectors.xy[MovingWindows.detectors.xy[ ,"upper.x",t] > xmax, "upper.x", t] <- xmax
     MovingWindows.detectors.xy[MovingWindows.detectors.xy[ ,"upper.y",t] > ymax, "upper.y", t] <- ymax
     
     MovingWindows.AC.xy[MovingWindows.AC.xy[ ,"lower.x",t] < xmin, "lower.x", t] <- xmin
     MovingWindows.AC.xy[MovingWindows.AC.xy[ ,"lower.y",t] < ymin, "lower.y", t] <- ymin
     MovingWindows.AC.xy[MovingWindows.AC.xy[ ,"upper.x",t] > xmax, "upper.x", t] <- xmax
     MovingWindows.AC.xy[MovingWindows.AC.xy[ ,"upper.y",t] > ymax, "upper.y", t] <- ymax
     
  }
  
  
  xy.bounds.detectors <- array(NA, c(n.individuals, 2, 2, n.years))
  xy.bounds.detectors[ , ,1, ] <- MovingWindows.detectors.xy[,c("lower.x", "lower.y"),]
  xy.bounds.detectors[ , ,2, ] <- MovingWindows.detectors.xy[,c("upper.x", "upper.y"),]
  
  
  xy.bounds.AC <- array(NA, c(n.individuals, 2, 2, n.years))
  xy.bounds.AC[ , ,1, ] <- MovingWindows.AC.xy[,c("lower.x", "lower.y"),]
  xy.bounds.AC[ , ,2, ] <- MovingWindows.AC.xy[,c("upper.x", "upper.y"),]
  
  
  # ---- STEP 5: IDENTIFY DETECTORS WITHIN EACH MOVING WINDOW ----- 
  detector.index.list <- list()
  n.detectors.id <- matrix(NA, n.individuals, n.years)
  for(t in 1:n.years){	
     detector.index.list[[t]] <- lapply(1:n.individuals, function(x){which(    detector.xy[ ,1,t] >= MovingWindows.detectors.xy[x,3,t]
                                                                               & detector.xy[ ,2,t] >= MovingWindows.detectors.xy[x,4,t]
                                                                               & detector.xy[ ,1,t] <= MovingWindows.detectors.xy[x,5,t]
                                                                               & detector.xy[ ,2,t] <= MovingWindows.detectors.xy[x,6,t])})
     n.detectors.id[ ,t] <- unlist(lapply(detector.index.list[[t]], length))
  }#t
  
  n.detectors.MW <- max(n.detectors.id)
  detector.index <- array(NA,c(n.individuals, n.detectors.MW, n.years))
  
  for(t in 1:n.years){
     detector.index[ , ,t] <- do.call(rbind, lapply(detector.index.list[[t]], function(x) c(x, rep(NA, n.detectors.MW - length(x)))))
  }#t
  
  
  # ---- STEP 6: KEEP ONLY MOVING WINDOWS WITH AT LEAST ONE DETECTOR AND ONE HABITAT ----- 
  # IF NO DETECTORS IN DETECTOR WINDOW  
  #non.null.det <- which(apply(detector.index, c(1), function(x)any(x>0)))
  non.null <- which(apply(detector.index, c(1), function(x)any(x>0)))
  
 

   # SUBSET THE DIFFERENT OBJECTS 
  n.individuals <- length(non.null)
  n.individuals.augmented <- n.individuals - n.individuals.detected
  
  detector.index <- detector.index[non.null, ,]
  n.detectors.id <- n.detectors.id[non.null,]

  MovingWindows.detectors.xy <- array(MovingWindows.detectors.xy[non.null, , ], c(n.individuals, 6, n.years))
  MovingWindows.AC.xy <- array(MovingWindows.AC.xy[non.null, , ], c(n.individuals, 6, n.years))
  
  xy.bounds.AC <- array(xy.bounds.AC[non.null, , , ], c(length(non.null), 2, 2, n.years))
  xy.bounds.detectors <- array(xy.bounds.detectors[non.null, , , ], c(length(non.null), 2, 2, n.years))
  
  
  dimnames(MovingWindows.AC.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  dimnames(MovingWindows.detectors.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  
  
  # ---- STEP 6: IDENTIFY DETECTIONS OUTSIDE EACH MOVING WINDOW ----- 
  outside.index.list <- id.outside <- dist.outside.id <- list()
  n.outside.id <- matrix(NA, n.individuals.detected, n.years)
  
  for(t in 1:n.years){
    id.detectors[[t]] <- lapply(1:n.individuals.detected, function(x){which(y[x, ,t] %in% observable.states)})
    outside.index.list[[t]] <- lapply(1:n.individuals.detected, function(x){id.detectors[[t]][[x]][which(!(id.detectors[[t]][[x]] %in% detector.index.list[[t]][[x]]))]})
    
    if(any(lapply(outside.index.list[[t]],length) > 0)){
      dist.outside.id[[t]] <- lapply(1:n.individuals.detected, function(i){
        do.call(c, lapply(outside.index.list[[t]][[i]],function(x){
          sqrt((detector.xy[x,1,t] - MovingWindows.detectors.xy[i,1,t])^2 + (detector.xy[x,2,t] - MovingWindows.detectors.xy[i,2,t])^2)
        }))
      })
    }#if
    
    n.outside.id[ ,t] <- unlist(lapply(outside.index.list[[t]], length))
    id.outside[[t]] <- which(n.outside.id > 0)
    
    if(length(id.outside[[t]]) > 0){
      writeLines(c( "WARNING :", "Some detections fall outside of the individual Moving Window",""))
      
      lapply(id.outside[[t]], function(x){
        temp <- do.call(paste, as.list(outside.index.list[[t]][[x]]))
        temp1 <- do.call(paste, as.list(round(dist.outside.id[[t]][[x]],2)))
        
        writeLines(c( paste(" The DETECTORS :  ", temp," are outside the Moving Window for INDIVIDUAL :", x)
                      , paste(" respectively at a DISTANCE of : ", temp1,"units from the center of the Moving Window")
                      , "")) 
      })
    }#if
  }#t
  
   # ---- STEP 7: OBTAIN PERCENT OF HABITAT IN EACH MOVING WINDOW ----- 
  
  dim.mx <- xy.bounds.AC
  dim.mx[,1,1,1] <-  dim.mx[,1,1,1] +1
  dim.mx[,2,1,1] <-  dim.mx[,2,1,1] +1
  
  # Obtain number of cells covered by window size 
  size.habitat.window <-  abs(x.seq[1]-x.seq[2])
  
  # Obtain percent of habitat cells covered by window size 
  prop.habitat.window <- rep(0,dim(dim.mx)[1])
  n.habitat.cells <- rep(0,dim(dim.mx)[1])
  prop.cells <- rep(0,dim(dim.mx)[1])
  prop.cells.searched <- rep(0,dim(dim.mx)[1])
  
  Cells.In.AC.Window <- matrix(NA, nrow=dim(dim.mx)[1], ncol=(size.habitat.window^2) )
  Cells.Searched.In.AC.Window <- matrix(NA, nrow=dim(dim.mx)[1], ncol=(size.habitat.window^2) )
  
  
   for(i in 1:dim(dim.mx)[1]){
    # define X an Y extent of the moving window 
    min.x <- dim.mx[i,1,1,1]
    max.x <- dim.mx[i,1,2,1]
    min.y <- dim.mx[i,2,1,1]
    max.y <- dim.mx[i,2,2,1]
    
    window.mx.ind <- as.matrix(habitat.mx[c(min.y:max.y), c(min.x:max.x)])
    window.ID.ind <- as.matrix(IDCells.mx[c(min.y:max.y), c(min.x:max.x)])
    
    if(sum(window.mx.ind)==0){
       next
    }
    
    if(augment.searched == TRUE){
    window.search.ind <- as.matrix(searched.mx[c(min.y:max.y), c(min.x:max.x)])
    Cells.Searched.In.AC.Window[i,1:length(window.ID.ind[window.search.ind==1])] <- window.ID.ind[window.search.ind==1]  
    prop.cells.searched[i] <- sum(!is.na(window.ID.ind[window.search.ind==1] ))/ (size.habitat.window^2)
    }
    
    Cells.In.AC.Window[i,1:length(window.ID.ind[window.mx.ind==1])] <- window.ID.ind[window.mx.ind==1]
    
    prop.habitat.window[i] <- length(window.ID.ind[window.mx.ind==1])/(size.habitat.window^2)
    n.habitat.cells[i] <- length(window.ID.ind[window.mx.ind==1])
    prop.cells[i] <- dim(window.ID.ind)[1]*dim(window.ID.ind)[2]/(size.habitat.window^2)
  }
  
  dimnames(MovingWindows.AC.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  dimnames(MovingWindows.detectors.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  
 
  # ---- STEP 8: SUBSET TO WINDOWS WITHIN SUITABLE HABITAT ----- 
  id.hab <- which(prop.habitat.window!=0)
  
  
  n.individuals <- length(id.hab)
  n.individuals.augmented <- n.individuals - n.individuals.detected
  
  detector.index <- detector.index[id.hab,]
  n.detectors.id <- n.detectors.id[id.hab]
  
  MovingWindows.detectors.xy <- array(MovingWindows.detectors.xy[id.hab, , ], c(n.individuals, 6, n.years))
  MovingWindows.AC.xy <- array(MovingWindows.AC.xy[id.hab, , ], c(n.individuals, 6, n.years))
  
  xy.bounds.AC <- array(xy.bounds.AC[id.hab, , , ], c(length(id.hab), 2, 2, n.years))
  xy.bounds.detectors <- array(xy.bounds.detectors[id.hab, , , ], c(length(id.hab), 2, 2, n.years))
  
  prop.habitat.window <- prop.habitat.window[id.hab]
  n.habitat.cells <- prop.habitat.window[id.hab]
  prop.cells <- prop.habitat.window[id.hab]
  prop.cells.searched <- prop.cells.searched[id.hab]
  
  Cells.In.AC.Window <- Cells.In.AC.Window[id.hab,]
  n.Cells.In.AC.Window <- apply(Cells.In.AC.Window, 1, function(x) sum(!is.na(x)))
  
  dimnames(MovingWindows.AC.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  dimnames(MovingWindows.detectors.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
  
  #Replace all the NAS by this value, so jags doesnt try to replace the NAs.
  Cells.In.AC.Window[is.na(Cells.In.AC.Window)] <- -999
  
  
  
  #AUGMENT MORE IN AREAS NOT SEARCHED 
  if(augment.searched==TRUE){
   nlayers.augm <- round(n.layers1 + prop.cells.searched * - (n.layers1-min.augment.layer))
   #nlayers.augm <- round(n.layers1 + exp(1-prop.cells.searched) * min.augment.layer)
   nlayers.augm[1:n.individuals.detected] <- 1
   prop.cells.searched1 <- prop.cells.searched
   plot(nlayers.augm~prop.cells.searched1, xlab="Prop cells searched", ylab="number of augmented layers")
   
   augment.id <- rep(1:n.individuals, nlayers.augm)
   
   n.individuals <- length(augment.id)
   n.individuals.augmented <- n.individuals - n.individuals.detected
   
   detector.index <- detector.index[augment.id,]
   n.detectors.id <- n.detectors.id[augment.id]
   
   MovingWindows.detectors.xy <- array(MovingWindows.detectors.xy[augment.id, , ], c(n.individuals, 6, n.years))
   MovingWindows.AC.xy <- array(MovingWindows.AC.xy[augment.id, , ], c(n.individuals, 6, n.years))
   
   xy.bounds.AC <- array(xy.bounds.AC[augment.id, , , ], c(length(augment.id), 2, 2, n.years))
   xy.bounds.detectors <- array(xy.bounds.detectors[augment.id, , , ], c(length(augment.id), 2, 2, n.years))
   
   prop.habitat.window <- prop.habitat.window[augment.id]
   n.habitat.cells <- prop.habitat.window[augment.id]
   prop.cells <- prop.habitat.window[augment.id]
   prop.cells.searched <- prop.cells.searched[augment.id]
   
   Cells.In.AC.Window <- Cells.In.AC.Window[augment.id,]
   n.Cells.In.AC.Window <- apply(Cells.In.AC.Window, 1, function(x) sum(!is.na(x)))
   
   dimnames(MovingWindows.AC.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
   dimnames(MovingWindows.detectors.xy) <- list(1:n.individuals, c("x", "y", "lower.x", "lower.y", "upper.x", "upper.y"), 1:n.years)
   }
  # ---- STEP 10: PLOT ----- 
  if(plot.check){
    ## just see an example for 10 individuals randomly selected for window size AC and window size detectors 
    t=1
    habitat.sp <- data.frame(habitat.xy[,,t]) 
    habitat.sp[ ,2] <- -habitat.sp[ ,2]
    coordinates(habitat.sp) <- habitat.sp 
    detector.sp <- data.frame(detector.xy[,,t]) 
    detector.sp[ ,2] <- -detector.sp[ ,2]
    coordinates(detector.sp) <- detector.sp
    
    par(mfrow=c(1,3), mar=c(5,1,10,1))
    
    plot(habitat.sp, col = "pink", pch = 19, cex = 0.1, main= "AC and detector windows for 10 \n randomly sampled detected individuals" )
    col <- rainbow(n.individuals.detected)
    col <- sample(col)
    plot(detector.sp, col =adjustcolor("blue",alpha.f = 0.2) , pch = 19, cex = 0.1, add = TRUE)
    
    lapply( sample(n.individuals.detected,10), function(i){                    
      xv <- c(MovingWindows.AC.xy[i,"upper.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"upper.x",t])
      yv <- c(-MovingWindows.AC.xy[i,"upper.y",t], -MovingWindows.AC.xy[i,"upper.y",t], -MovingWindows.AC.xy[i,"lower.y",t],-MovingWindows.AC.xy[i,"lower.y",t])
      polygon(xv, yv, border = "gray80", col = adjustcolor(col[i],alpha.f = 0.8))
      
      xv.det <- c(MovingWindows.detectors.xy[i,"upper.x",t], MovingWindows.detectors.xy[i,"lower.x",t], MovingWindows.detectors.xy[i,"lower.x",t], MovingWindows.detectors.xy[i,"upper.x",t])
      yv.det <- c(-MovingWindows.detectors.xy[i,"upper.y",t], -MovingWindows.detectors.xy[i,"upper.y",t], -MovingWindows.detectors.xy[i,"lower.y",t],-MovingWindows.detectors.xy[i,"lower.y",t])
      polygon(xv.det, yv.det, border = "gray80", col = adjustcolor(col[i],alpha.f = 0.4))
      
    })
    
    #for all individuals plot AC regions 
    for(t in 1:n.years){
      
      habitat.sp <- data.frame(habitat.xy[,,t]) 
      habitat.sp[ ,2] <- -habitat.sp[ ,2]
      coordinates(habitat.sp) <- habitat.sp 
      detector.sp <- data.frame(detector.xy[,,t]) 
      detector.sp[ ,2] <- -detector.sp[ ,2]
      coordinates(detector.sp) <- detector.sp
      
      # ----- DETECTED
      plot(habitat.sp, col = "pink", pch = 19, cex = 0.1, main= "AC window for detected individuals")
      lapply(1:n.individuals.detected, function(i){                    
        xv <- c(MovingWindows.AC.xy[i,"upper.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"upper.x",t])
        yv <- c(-MovingWindows.AC.xy[i,"upper.y",t], -MovingWindows.AC.xy[i,"upper.y",t], -MovingWindows.AC.xy[i,"lower.y",t],-MovingWindows.AC.xy[i,"lower.y",t])
        polygon(xv, yv, border = "gray80", col = grey(0.5, alpha = 0.5))
      })
      plot(detector.sp, col = "blue", pch = 19, cex = 0.1, add = TRUE)
      
      # ----- OUTSIDE
      lapply(id.outside[[t]], function(i){                    
        xv <- c(MovingWindows.detectors.xy[i,"upper.x",t], MovingWindows.detectors.xy[i,"lower.x",t], MovingWindows.detectors.xy[i,"lower.x",t], MovingWindows.detectors.xy[i,"upper.x",t])
        yv <- c(-MovingWindows.detectors.xy[i,"upper.y",t], -MovingWindows.detectors.xy[i,"upper.y",t], -MovingWindows.detectors.xy[i,"lower.y",t], -MovingWindows.detectors.xy[i,"lower.y",t])
        polygon(xv, yv, border = "red", col = rgb(1, 0, 0, alpha = 0.1))
        plot(detector.sp[outside.index.list[[t]][[i]], ], col = "red", pch = 19, cex = 0.5, add = TRUE)
        lapply(outside.index.list[[t]][[i]], function(j){
          segments( x0 = detector.sp$X1[j], y0 = detector.sp$X2[j]
                    , x1 = MovingWindows.detectors.xy[i,"x",t], y1 = -MovingWindows.detectors.xy[i,"y",t],col = "red")})
      })
      
      # ----- AUGMENTED
      plot(habitat.sp, col = "pink", pch = 19, cex = 0.1, main= "AC window for augmented individuals")
      lapply((n.individuals.detected+1):n.individuals, function(i){                    
        xv <- c(MovingWindows.AC.xy[i,"upper.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"lower.x",t], MovingWindows.AC.xy[i,"upper.x",t])
        yv <- c(0-MovingWindows.AC.xy[i,"upper.y",t], 0-MovingWindows.AC.xy[i,"upper.y",t], 0-MovingWindows.AC.xy[i,"lower.y",t], 0-MovingWindows.AC.xy[i,"lower.y",t])
        polygon(xv, yv, border = "white", col = grey(0.5, alpha = 0.3))
      })
      plot(detector.sp, col = "blue", pch = 19, cex = 0.1, add = TRUE)
    }#t
    
    if(augment.searched==TRUE){
       par(mar=c(4,4,4,4))
       plot(nlayers.augm~prop.cells.searched1, xlab="Prop cells searched", ylab="number of augmented layers")
    }
    
    
    
    
  }#if
  
  # ---- STEP 11: OUTPUT -----
  out <- list(   y.augmented = abind(y, array(0, c(n.individuals.augmented, n.detectors, n.years)), along = 1)
               , n.detectors = n.detectors.id
               , detector.index = detector.index
               , xy.bounds.AC = xy.bounds.AC
               , xy.bounds.detectors = xy.bounds.detectors
               , MovingWindows.detectors.xy = MovingWindows.detectors.xy
               , MovingWindows.AC.xy = MovingWindows.AC.xy
               , n.individuals = n.individuals
               , n.individuals.detected = n.individuals.detected
               , n.individuals.augmented = n.individuals.augmented
               , augmented.id.density = n.layers/(window.AC*window.AC)
               , prop.habitat.window = prop.habitat.window
               , n.habitat.cells=n.habitat.cells
               , Cells.In.AC.Window = Cells.In.AC.Window
               , n.Cells.In.AC.Window = n.Cells.In.AC.Window
               , prop.cells=prop.cells
               , number.detections.outside.window = n.outside.id)
  
  return(out)
}
