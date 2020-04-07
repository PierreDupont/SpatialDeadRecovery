#' @title GetPopCompo extracts the state-composition of the population from the z state-matrix
#'
#' @description
#' \code{GetPopCompo} returns a  \code{list} object with the composition of the population (number of individuals in each state) at each occasion
#'
#' @param z A \code{matrix} of dimensions number of individuals*number of years containing the individual states 
#'	@param plot.check A \code{logical} to display validation plots.
#'	
#' @return A \code{list} object with the composition of the population (number of individuals in each state) at each occasion
#' 
#' @examples
#' z <- matrix(sample(5, 500, replace=TRUE), 50, 10)
#' myPopCompo <- GetPopCompo(z = z)

GetPopCompo <- function( z
                       , plot.check = TRUE
                       , cex = 1.5
                       , lwd = 2
                       ,legend.placement="topright")
{
   ## IDENTIFY THE DIFFERENT STATES IN z
   z.levels <- unique(na.omit(unlist(apply(z, 1, unique))))  
   z.levels <- z.levels[order(z.levels)]
   
   ## COUNT THE NUMBER OF INDIVIDUALS PER STATE (except NAs)
   Pop.Compo <- list()
   for(l in 1:length(z.levels)){
      Pop.Compo[[l]] <- apply(z, 2, function(x){length(which(x == z.levels[l]))})
      }#l
   Pop.Compo[[length(z.levels) + 1]] <- apply(z, 2, function(x){length(na.omit(x))})
   
   ## PLOT POPULATION DYNAMICS
   if(plot.check){
      plot(1:dim(z)[2], Pop.Compo[[1]], type="l", col=1, lwd=2, ylim=c(0,max(Pop.Compo[[length(z.levels)+1]])), xlab = "Years", ylab = "# of Individuals")
      for(l in 2:length(Pop.Compo)){
         points(1:dim(z)[2],Pop.Compo[[l]], type="l", col=l, lwd=2)
         }#l
      legend(legend.placement
            , legend =  c(unlist(lapply(z.levels, function(x){paste(" State", x)})), "Total")
            , col = 1:length(Pop.Compo)
            , lwd = lwd
            , cex = cex
            , inset = 0.01)
      }#if
   
   return(Pop.Compo)
}


