#' @title GetPopDensity is a function to extract the individual and population AC densities.
#'
#' @description
#' \code{GetPopDensity} returns a \code{list} object with diferent rasters of individuals and population's density
#' 
#' @param habitat.mx A \code{matrix} object containing the id of the habitat cells
#'	@param r.origin A \code{raster} object from makeHabitat function 
#' @param posterior.sx A \code{array} object from the jags output containing the x-coordinates of the individual ACs 
#' @param posterior.sy A \code{array} object from the jags output containing the y-coordinates of the individual ACs 
#' @param posterior.z A \code{array} object from the jags output containing the individual states 
#' @param alive.states A \code{vector} indicating which states in the posterior.z object correspond to alive individuals
#' @param time.interval A \code{vector} with the indices of the years to be kept
#' @return A \code{list} with :
#' PopCounts.ar: an array with the total number of ACs per habitat cell (dimensions: 1.y-coordinates; 2.x-coordinates; 3.years)
#' PopDens.ar: an array with the mean number of Acs per habitat cell
#' PopCounts.r.stack: a raster stack with the total number of ACs per habitat cell
#' PopDens.r.stack: a raster stack with the mean number of Acs per habitat cell
#' 
#' 
#' @examples
#' myPopDensity <- GetDensityMetrics( habitat.mx = myHabitat.list$habitat.mx
#'                                  , r.origin= myHabitat.list$habitat.r
#'                                  , posterior.sy = my.jags.output$sim.list$sxy[,,2]
#'                                  , posterior.sx = my.jags.output$sim.list$sxy[,,1]
#'                                  , posterior.z = my.jags.output$sim.list$z
#'                                  , alive.states = 1
#'                                  , time.interval = NULL)


GetPopDensity <- function( habitat.mx                 # 0/1 matrix of the habitat
                         , r.origin                   # Original raster of the habitat 
                         , posterior.sy               # Posterior y coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
                         , posterior.sx               # Posterior x coordinate of individual ACS returned by jags (three dimensional: 1=iterations; 2=individuals; 3=time)
                         , posterior.z                # Posteriors state z returned by jags (three dimensions: 1=iterations; 2=individuals; 3=time)
                         , alive.states = NULL        # Alive states in z
                         , time.interval = NULL)      # Index of years to be returned
{
   ## Set the parameters of the results to be returned
   y.max <- dim(habitat.mx)[1]
   x.max <- dim(habitat.mx)[2]
   n.iter <- dim(posterior.z)[1]
   n.individuals <- dim(posterior.z)[2]
   n.years <- ifelse(length(dim(posterior.z))==2, 1, ifelse(is.null(time.interval),dim(posterior.z)[3],length(time.interval)))
   if(is.null(time.interval)){time.interval <- 1:n.years}
   
   ## Filter the posteriors
   sx <- sy <- Z <- array(0, dim = c(n.iter, n.individuals, n.years))
   if(length(dim(posterior.sy))!=2){sy <- posterior.sy[ , ,time.interval]
                                    sx <- posterior.sx[ , ,time.interval]}
   if(length(dim(posterior.sy))==2){
      for(t in 1:n.years){
         sy[,,t] <- posterior.sy
         sx[,,t] <- posterior.sx
         }
      }
   
   ## Assign a 1 to all alive individuals
   Z[which(posterior.z %in% alive.states)] <- 1
   
   # Transform AC locations to cell ID
   sx <- trunc(sx)+1
   sy <- trunc(sy)+1
   
   ## Keep AC locations, Iterations and time of alive individuals only (z == 1)
   SY <- sy[Z==1]
   SX <- sx[Z==1]
   TIME <- which(Z == 1, arr.ind = TRUE)[ ,3] 
   
   ## Create arrays of Population counts and density 
   
   #x<-3
   AC.Pop.Counts <- array(0, c(y.max, x.max, n.years))
   for (x in 1:length(SX)){
      AC.Pop.Counts[SY[x],SX[x],TIME[x]] <- AC.Pop.Counts[SY[x],SX[x],TIME[x]] + 1
      }
   AC.Pop.Density <- AC.Pop.Counts/n.iter
   
   ## Create raster stacks of yearly population densities and individuals counts
   for(t in 1:n.years){
      r.Pop.Density <- raster(AC.Pop.Density[,,t])
      r.Pop.Counts <- raster(AC.Pop.Counts[,,t])
      extent(r.Pop.Density) <- extent(r.origin)
      extent(r.Pop.Counts) <- extent(r.origin)
      res(r.Pop.Density) <- res(r.origin)
      res(r.Pop.Counts) <- res(r.origin)
      if(t>1){PopDens.r.stack <- stack(PopDens.r.stack, r.Pop.Density)}else{PopDens.r.stack <- r.Pop.Density}
      if(t>1){PopCounts.r.stack <- stack(PopCounts.r.stack, r.Pop.Counts)}else{PopCounts.r.stack<- r.Pop.Counts}
      }
   
 ## Create Output list
 return(list( PopCounts.ar = AC.Pop.Counts
            , PopDens.ar = AC.Pop.Density
            , PopCounts.r.stack = PopCounts.r.stack
            , PopDens.r.stack = PopDens.r.stack))
   }
