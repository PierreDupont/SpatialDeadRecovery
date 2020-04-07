#' CECI EST UNE DESCRIPTION
#' @title Data set simulation function.
#'
#' @description
#' \code{SimulateZ.PopDyn} returns a \code{list} object with the composition of the simulated population,
#'  its pop size and the corresponding z matrix for CMR models.
#'
#' @param N0 A \code{numeric} object containig the original number of individuals in the pop to be simulated (alternative to z0)
#'	@param n.occasions A \code{numeric} object containing the number of years tobe simulated.
#' @param PHI A \code{matrix} object containing the state-specific survival probabilities.
#' @param REPRO A \code{matrix} object containing the state-specific reproduction probabilities.
#' @param FEC A \code{vector} object containing the state-specific fecundity (= mean litter size).
#' @param INIT.STATE A \code{vector} object containing the initial-state probabilities (= prob to enter the data set in each state).
#' @param z0 A \code{vector} object containing the state of the different individuals at t0 (alternative to N0).

#' @return A \code{list} object with the population attributes 
#' pop.ls: a list of population composition for each time step
#' N: a list of pop sizes
#' z: the state matrix to be used in CMR models
#' 
#' @examples
#' Generate simulated population with three states: 1=JUVENILES, 2=ADULTS, 3=DEAD
#' PHI <- matrix(c(0,0.85,0.15,
#'                 0,0.75,0.25,
#'                 0,0,1),n.row=3,byrow=TRUE)
#' REPRO <- matrix(c(0,0,1,
#'                   0,0.75,0.25,
#'                   0,0,1),n.row=3,byrow=TRUE)
#' FEC <- c(0,1,0)
#' mySim <- SimulateZ.PopDyn( N0 = NULL
#'                          , n.occasions = 5
#'                          , PHI = PHI
#'                          , REPRO = REPRO
#'                          , FEC = FEC
#'                          , INIT.STATE = NULL
#'                          , z0 = c(1,1,1,1,2,1,1,2,2,2)

SimulateZ.PopDyn <- function( N0 = NULL
                            , n.occasions
                            , PHI
                            , REPRO
                            , FEC  
                            , INIT.STATE = NULL
                            , z0 = NULL)
{
   ## ---- CREATE POPULATION LIST ----
   if(is.null(N0)){N0 <- length(z0)}
   if(is.null(z0)){z0 <- sample(dim(PHI)[1],N0,replace=TRUE)}
   POP <- POP.SIZE <- list()
   POP[[1]] <- data.frame(id=1:N0, z=z0)
   
   ## ---- CREATE POPULATION SIZE VECTOR ----
   N <- vector()
   N[1] <- N0
   
   #-----------------------------------------------------------------------
   ## ---- LOOP POPULATION PROCESS OVER n.occasions ----
   for(t in 2:(n.occasions+1))
      {
      Z <- R <- B <- Z.new <-  NULL
      for (i in 1:(dim(POP[[t-1]])[1]))
         {
         Z[i] <- which(rmultinom(1, 1, PHI[POP[[t-1]]$z[i], ])==1)
         R[i] <- which(rmultinom(1, 1, REPRO[Z[i], ])==1)
         B[i] <- rpois(1, FEC[R[i]])
         }
      N.new <- sum(B[])
      if(N.new>=1)
         {
         if(is.null(INIT.STATE)){INIT.STATE <- c(1,rep(0,dim(PHI)[1]-1))}
         for (i in 1:N.new)
            {
            Z.new[i] <- which(rmultinom(1, 1, INIT.STATE)==1)
            }
         Z <- c(Z,Z.new)
         }
      N[t] <- length(Z)
      POP[[t]] <- data.frame(id=1:N[t], z = Z)
      }
   #-----------------------------------------------------------------------
   ## ---- GENERATE z BASED ON POPULATION LIST
   z <- matrix(NA, N[n.occasions+1], n.occasions+1)
   for (t in 1:(n.occasions+1))
      {
      z[1:N[t],t] <- POP[[t]]$z
      pop.size <- vector()
      for (s in 1:dim(PHI)[1])
         {
         pop.size[s] <- length(which(z[ ,t]==s))
         }
      POP.SIZE[[t]] <- c(pop.size, N[t])
      }
   
   #----------------------------------------------------------------------
   return (list(pop.ls = POP, N = POP.SIZE, z = z))
   }