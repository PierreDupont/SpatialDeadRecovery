#' @title Function to reduce a spatial lines object to regularly spaced points (based on a grid).
#'
#' @description
#' \code{Tracks2Points} returns a SpatialPointsDataframe object with \code{}. The data attribute of the output file contains a field specifying the total length of segments that are depicted by a point at that location.
#' 
#' @param track SpatialLinesDataframe.
#' @param poly Polygon for spatial subsetting of \code{track}.
#' @param resolution Resolution of the discretized output (in distance units of \code{track}).
#' @param plot Logical for whether (\code{TRUE}) or not (\code{FALSE}) plots are to be generated.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/Tracks2Points.R
#' @keywords simul
#'
#' @examples
#' # Convert a spatial lines dataframe to regular spatial points dataframe with custom spacing:
#' track.sp<-Tracks2Points(track=trail.sl,poly=area.poly,resolution=250,plot.check=TRUE)


Tracks2Points <- function(track, 
                          polygon = NULL, 
                          resolution, 
                          plot.check = TRUE)
   {
   projection(track) <- projection(polygon)                  #---just in case...
   
   if(!is.null(polygon)){track <- intersect(track,polygon)}  #---clip to desired polygon, if provided
   
   r <- raster(extent(track), resolution=resolution)
   projection(r) <- projection(track)                        #---just in case...
   r <- rasterize(x=track, y=r,fun='length')
   
   track.xy <- data.frame(xyFromCell(r, 1:ncell(r)))
   track.xy$length <- r[]
   track.xy <- na.omit(track.xy)
   
   track.sp <- SpatialPointsDataFrame(track.xy[,c("x","y")], data=track.xy, proj4string=CRS(projection(track)))
   
   if(plot.check == TRUE)
      {
      plot(r)
      plot(track, add=TRUE)
      plot(track.sp, add=TRUE, col="red", pch=19, cex=0.6)
      if(!is.null(polygon))plot(polygon,add=TRUE,lwd=2)
      }
  
   return(list(track.sp = track.sp, track.r = r, track.lines = track))
}
