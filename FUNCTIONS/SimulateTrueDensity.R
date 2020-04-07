#' @title Function to simulate the true density from the habitat raster from half normal distribution
#'
#' @description
#' \code{SimulateTrueDensity} simulates detection of individuals (acitivity centers) from a Poisson or Bernoulli(Binomial) distribution.
#'  it returns the true density in the study area given some simulations activity centers and a sigma and lambda/p0
#'  it also has the possibility to repeat the procedure N number of times in order to approach the true density.
#'    
#' @param myHabitat.list a list of habitat file from MakeHabitat function 
#' @param lambda0 Numeric variable denoting the intercept value of the negative binomial detection function describing the decay of detection with increasing distance from the AC.
#' @param sigma Numeric variable denoting the scale parameter of the detection function.
#' @param AC.sp Spatial points dataframe with individual activity center locations.
#' @param type String the type of distribution to be used c("Poisson", "Binomial", "Bernoulli"). "Binomial" and "Bernoulli" are the same. 
#' @param rep Numeric for the number of times the detection of individuals should be applied. 
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 
#' @keywords simul
#'
#' @examples
#'grid.size <- 32
#'coords = matrix(c(0        , 0 ,
#'                  grid.size, 0 ,
#'                  grid.size, grid.size,
#'                  0        , grid.size,
#'                  0        ,  0 
#'), ncol = 2, byrow = TRUE)
#'
#'P1 <- Polygon(coords)
#'myStudyArea.poly <-  SpatialPolygons(list(Polygons(list(P1), ID = "a")),
#'                                     proj4string=CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#'
#'## Generate habitat matrix 
#'myHabitat.list <- MakeHabitat( poly = myStudyArea.poly                      ## Polygon of the focal study area
#'                               , resolution = 1                               ## Grid cell size (in meters)
#'                               , buffer = 0                               ## Size of the buffer area (in meters)
#'                               , polygon.clip = NULL                        ## Polygon of suitable habitat 
#'                               , plot.check = TRUE)   
#'
#' Simulate AC centers 
#'AC.sp <- spsample( myHabitat.list$buffered.habitat.poly, 20 ,type="random")
#'AC.sp$ID <- 1:length(AC.sp)
#'plot(myHabitat.list$buffered.habitat.poly)
#'points(AC.sp, pch=16, col="red")
#'  
#' true.density <-  SimulateTrueDensity(myHabitat.list = myHabitat.list
#'                                      , AC.sp = AC.sp
#'                                      , lambda0 = 10
#'                                      , sigma = 3
#'                                      , type = "Poisson"
#'                                      , rep = 50
#'                                      , plot.check = TRUE)
#' sum(values(true.density))# should be equal to the number of id simulated # 
SimulateTrueDensity <- function( myHabitat.list
                               , AC.sp
                               , sigma
                               , plot.check = TRUE
                               , study.area = NULL
                               , breaks = NULL)
   {
   ## Compute the distances between each habitat cells and individual ACs ====
   row.names(AC.sp) <- 1:length(AC.sp)
   D <- gDistance(myHabitat.list$habitat.sp, AC.sp, byid = TRUE)
   
   ## Calculate space-use (i.e. percentage of time spent in a given habitat cell)
   lambda <- exp(-D*D/(2*sigma*sigma))

   ## Obtain the cells that are unsuitable habitat
   cells <- which(values(myHabitat.list$habitat.r)==0) 
   
   ## Convert space-use in a raster for the first id
   r <- raster(matrix(lambda[1,], ncol=ncol(myHabitat.list$habitat.r), nrow=nrow(myHabitat.list$habitat.r), byrow = TRUE))
   extent(r) <- extent(myHabitat.list$habitat.r)
   res(r) <- myHabitat.list$resolution
   rval <- values(r)
   rval[cells] <- NA
   values(r) <- rval/ sum(rval, na.rm=T)
   
   ## Convert space-use in a raster for all other ids  
   for(i in 2:length(AC.sp)){
      r1 <- raster(matrix(lambda[i,], ncol=ncol(myHabitat.list$habitat.r), nrow=nrow(myHabitat.list$habitat.r), byrow = TRUE))
      extent(r1) <- extent(myHabitat.list$habitat.r)
      res(r1) <- myHabitat.list$resolution
      r1val <- values(r1)
      r1val[cells] <- NA
      values(r1) <- r1val/sum(r1val, na.rm = TRUE)
      r <- stack(r, r1)
      }
   
   ## Sum the raster stack to get total density 
   rtot <- sum(r)
   
   ## Plot it 
   if(plot.check==TRUE){
      if(is.null(breaks)){breaks <- seq(0, max(na.omit(values(rtot))), length.out = 20)} 
      col.pal <- rev(terrain.colors(length(breaks)))
      plot(rtot, breaks = breaks, col = col.pal)
      points(AC.sp, col="red", pch=16)
      if(!is.null(study.area)){plot(study.area, add=TRUE)}
      }
   
   return(rtot)
   } 



