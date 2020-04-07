#' @title Simulate a Markov Chain (e.g. to obtain a simulated state history).
#'
#' @description
#' \code{SimulateMC} returns a list object with \code{} 
#' 
#' @param t Numeric variable denoting the number of time steps.
#' @param P Matrix with transition probabilities.
#' @param x1 A numeric vector indicating the initial states.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SimulateMC.R
#' @keywords simul
#'
#' @examples
#' # Simulate an individual state history
#' trans.mx<-rbind(c(0.9,0.1,0),c(0,0.9,0.1),c(0,0,1))
#'init.state<-1
#'n.timesteps=100
#'SimulateMC(n.timesteps,trans.mx,init.state)
#' 



SimulateMC <- function(n.timesteps, trans.mx, init.state) {
   sim <- as.numeric(n.timesteps)
   m <- ncol(trans.mx)
   if (missing(init.state)) {
      sim[1] <- sample(1 : m, 1)
   } else {
      sim[1] <-  init.state
   }
   for (i in 2 : n.timesteps) {
      newstate <- sample(1 : m, 1, prob = trans.mx[sim[i - 1], ])
      sim[i] <- newstate
   }
   return(sim)
}

