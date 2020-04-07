#' @title GetTransitionSurfaces
#' 
#' @description
#' \code{GetTransitionSurfaces} Estimates surfaces of probability of success and its uncertainty (SD) for transition between states of departures and states of arrival from posterior z 
#'

#' @param habitat.mx A \code{matrix} of habitat (0/1) from MakeHabitat function
#' @param habitat.r A \code{raster} of habitat at the original scale from MakeHabitat function
#' @param IDCells.mx A \code{matrix} with the cell ID of the habitat MakeHabitat function 
#' @param posterior.sy A \code{matrix} of posterior sy of  jags. is typically a matrix with individuals (i) in columns i.e. sy[,i] 
#' @param posterior.sx A \code{matrix} of posterior sx of jags. is typically a matrix with individuals (i) in columns i.e. sx[,i] 
#' @param posterior.z A \code{matrix} of posteriors z of  jags. is typically a matrix with individuals (i) in columns i.e. z[,i] 
#' @param states.from A \code{vector} of states representing the states at the departure of the transition from which a probability of success of the transition will be calculated.  
#' @param states.to A \code{vector} of states representing the states at the arrival of the transition from which a probability of success of the transition will be calculated
#' @param cell.at.departure A \code{logical}. If TRUE, the location of the AC at the state of departure is used (AC[t]; states.from[t], states.to[t+1]); If FALSE, the location of the AC at the state of arrival is used (AC[t+1]; states.from[t], states.to[t+1]).
#' @param plot.check A \code{logical} if average surfaces over years should be plotted
#' 
#'
#' @return Transition surfaces rasters (mean, sd) for every transitions(stack) and averaged over all transitions (occasions)
#' @usage 
#' my.Transition.Surfaces <- GetTransitionSurfaces( habitat.mx  = myHabitat$habitat.mx               # 0/1 matrix of the study area
#'                                                , habitat.r   = myHabitat$habitat.r                # Original raster of the habitat 
#'                                                , IDCells.mx = myHabitat$IDCells.mx
#'                                                , posterior.sy = my.jags.output$JAGS.output$sims.list$sxy[,,2]            # Posterior y coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
#'                                                , posterior.sx = my.jags.output$JAGS.output$sims.list$sxy[,,1]            # Posterior x coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
#'                                                , posterior.z  =  my.jags.output$JAGS.output$sims.list$z          # Posteriors state z returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
#'                                                , states.from =  c(2) 
#'                                                , states.to =   c(3,4) 
#'                                                , cell.at.departure = T
#'                                                , plot.check = F)
#' @export
#' @examples
#' 
GetTransitionSurfaces <- function(
     habitat.mx                 # 0/1 matrix of the study area
   , habitat.r                   # Original raster of the habitat 
   , IDCells.mx  
   , posterior.sy               # Posterior y coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
   , posterior.sx               # Posterior x coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
   , posterior.z                # Posteriors state z returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
   , states.from =  c(2)        # states definition
   , states.to =   c(3,4)            
   , cell.at.departure = TRUE # if true, it means that transition surfaces are calculated after leaving the cell, if false, transiton calculted before the arrival to the cell.
   , plot.check = TRUE){
      
      
   #--- CHECK OBJECTS ##
    if(length(dim(posterior.z))<3){
         stop("Need more than one time step to compute transition surfaces!")                
    }
   
      
   #--- get some objects ready ##
  n.iter <- dim(posterior.z)[1]
  n.individuals <- dim(posterior.z)[2]
  n.years <- dim(posterior.z)[3]
  n.occasions <- n.years-1
  time.interval <- 1:n.years
  y.max <- dim(habitat.mx)[1]
  x.max <- dim(habitat.mx)[2]
  
  
  #--- OBJECT SET UP ##
  # IF THERE IS NO MOVEMENT IN ACS, ASSIGN SY AND SY TO EVERY YEAR
      sx <- array(0, dim = c(n.iter, n.individuals, n.years))
      sy <- array(0, dim = c(n.iter, n.individuals, n.years))
      
      if(length(dim(posterior.sy))!=2){
         sy <- posterior.sy
         sx <- posterior.sx
         }
      
      if(length(dim(posterior.sy))==2){
         for(t in 1:n.years){
            sy[,,t] <- posterior.sy
            sx[,,t] <- posterior.sx
         }
      }
      
   #--- ASSIGN POSTERIOR COORDINATES TO A X AND Y OF A CELL
   sx <- trunc(sx) + 1
   sy <- trunc(sy) + 1
   
#---- GET TRANSITION SURFACES ----
# create some temporary raster to store the output #  
mean.trans.r <- habitat.r
mean.trans.r[] <- NA
sd.trans.r <- mean.trans.r
rtmp <- mean.trans.r
rtmp1 <- mean.trans.r


if(cell.at.departure==T){
   occ <- c(1:(n.years-1))
   occ.t.plus1 <- c(2:n.years)
   
   statesfrom <- states.from  
   statesto <- states.to
   
}else{
   occ <- c(2:n.years)
   occ.t.plus1 <- c(1:(n.years-1))
   
   statesfrom <- states.to 
   statesto <- states.from 
   
}
   
   
# create an empty array to store the results 
arr <- array(NA, dim = c(y.max,x.max, n.iter, n.years-1))


   for(t in 1:n.years-1){# occasions 
      for(j in 1:n.iter){# iterations 
         
      # get which individual are in the states from at time t   
      id.alive.t <- which(posterior.z[j,,occ[t]] %in% statesfrom)
      
      # get their sxy coordinates 
      sx.alive.t <- sx[j, id.alive.t, occ[t]]
      sy.alive.t <- sy[j, id.alive.t, occ[t]]

      # get cell IDs of individual alive 
      cell <- numeric()
      for( id in 1:length(sx.alive.t)){
         cell[id] <-  IDCells.mx[sy.alive.t[id], sx.alive.t[id]]
      }
      
      # Check how many of individuals of each cell are in the state.to 
      cellID <- unique(cell)
      for(id in 1:length(cellID)){
         next.state <-  which(cell %in% cellID[id])
         arr[sy.alive.t[cell%in%cellID[id]][1], sx.alive.t[cell%in%cellID[id]][1], j, t] <- sum(posterior.z[j, id.alive.t[next.state], occ.t.plus1[t]] %in% statesto) / length(next.state)
      }#id
      #print(j)
      
            }#j

      rtmp[] <- apply(arr[,,,t] , 1:2, function(x) mean(x,na.rm=T))
      rtmp1[] <- apply(arr[,,,t] , 1:2, function(x) sd(x,na.rm=T))# try with adding a log transformation on the sd
      
      if(t>1){mean.trans.r <- stack(mean.trans.r,rtmp)}else{mean.trans.r[] <- rtmp[]}
      if(t>1){sd.trans.r <- stack(sd.trans.r, rtmp1)}else{sd.trans.r[] <- rtmp1[]}
      names(mean.trans.r) <- paste( rep("transition", length(names(mean.trans.r))), c(1:length(names(mean.trans.r))))
      names(sd.trans.r) <- paste( rep("transition", length(names(sd.trans.r))), c(1:length(names(sd.trans.r))))
   }

### get average transition over years # 
average.mean.trans.r <- mean(mean.trans.r ,na.rm=T)
average.sd.trans.r <- mean(sd.trans.r ,na.rm=T)

#---- PLOT CHECK ----
if(plot.check==TRUE){
par(mfrow=c(1,2))
plot(average.mean.trans.r, legend.args=list(text='transition probability'), axis.args=list(cex.axis=0.6)
     , main= paste("Average transition probability \n from state ", states.from, " to state ", states.to, sep="") )
plot(average.sd.trans.r, legend.args=list(text='transition probability'), axis.args=list(cex.axis=0.6)
     , main= paste("Average uncertainty (sd) in transition probability \n from state ", paste(states.from,sep=" , "), " to state ",
                   paste(states.to,sep=","), sep="") )
}


#---- RETURN SURFACES ----
retu <- list(  average.mean.trans.r=average.mean.trans.r
        , average.sd.trans.r=average.sd.trans.r
        , stack.mean.trans.r=mean.trans.r
        , stack.sd.trans.r = sd.trans.r) 

return(retu)


}



   
