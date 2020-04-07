#' @title Simulate a state matrix for multiple individuals over multiple years.
#'
#' @description
#' \code{SimulateMultiYearZ} returns a list object with \code{} 
#' 
#' @param n Numeric variable denoting the number of individuals.
#' @param ntime Numeric variable denoting the number of occasions.
#' @param trans Transition matrix.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/BuildTransitionMatrix.R
#' @keywords simul
#'
#' @examples
#' # Build an example transition matrix:
#' SimulateMultiYearZ(10,20,P)
#' 



SimulateMultiYearZ <- function(n, n.timesteps, trans.mx,init.state=rep(1,n)){
   ch <- matrix(NA, ncol = n.timesteps, nrow = n)
   for(i in 1:n){
      ch[i, ] <- SimulateMC(n.timesteps=n.timesteps, trans.mx=trans.mx, init.state =init.state[i])
   }
   return(ch)
   
}

