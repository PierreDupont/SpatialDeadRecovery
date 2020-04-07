#' @title State Matrix (z) Simulation function based on multi-state modelling processes
#'
#' @description
#' \code{SimulateZ.MultiState} returns a \code{matrix} object with the individual states for each time step.
#'
#' @param N0 A \code{vector} object containing the number of individuals to enter the data set at each time occasion.
#' @param n.individuals A \code{numeric} object to define the total number of individuals.
#' @param n.occasions A \code{numeric} object to define the number of capture occasions.
#'	@param init.state A \code{vector}, {matrix} or {array} object containing the probabilities to enter the data set in a given state. If covariates effect are added, state.trans should be a character matrix with all parameters added as string. i.e.("phi") 
#' @param state.trans A \code{matrix} or {array} object representing the state-transition matrix. If covariates are used, matrix should be filled with characters (see examples)
#' @param plot.check A \code{logical} to display z image and covariates effect
#' @param id.covs A \code{list} of \code{dataframe} with individual covariates for each  transition parameters. The name of each element of the list should match the names entered in the state.trans object
#' @param id.betas A \code{list} of \code{vector} with beta coefficients on the log scale for each parameter. The interecept used is the one entered in the \code{intercepts}. The parameters should match the names entered in the state.trans object
#' @param intercepts  A \code{list} of \code{vector} with intercept values to be used for each parameter. If intercept values differs between years, provide vector of intercept values for each year and each transition parameters.
#' @return A \code{matrix} object with individual state-histories z.
#' 
#' @examples
#' # with no covariates# 
#' N0 <- c(10,5,20,8)
#' init.state <- c(0.8,0.2,0)
#' state.trans <- matrix(c(0.3,0.5,0.2,0,0.8,0.2,0,0,1),3,3,byrow = TRUE)
#' SimulateZ.MultiState(N0 = N0,init.state = init.state,state.trans = state.trans)  
#' 
#' # With Covariates (phi/psi) and time varying intercept on PHI
#' state.trans <- matrix(c( "1-psi", "psi"  , 0           , 0                  ,
#'                 0    , "phi" , "(1-phi)*r "  ,"(1-phi)*(1-r)"    , 
#'                 0    , 0    , 0           , 1                  ,
#'                 0    , 0    , 0           ,                  1), nrow = 4, byrow = TRUE)
#' 
#' 
#' id.covs.phi <- data.frame(Roads = rnorm(300 ,mean = 0, sd=1),Forest = rnorm(300,mean = 0, sd=1))
#' id.covs.psi <- data.frame(Roads = id.covs.phi[,1])
#' id.betas.phi <- c(-0.3,0.2) # beta coefficient 
#' id.betas.psi <- c( -0.3) # beta for the covariates 
#' z.mx <- SimulateZ.MultiState( N0 = c(300,0,0,0), init.state = c(1,0,0,0),  state.trans = m , plot.check = TRUE, id.covs = list( phi=id.covs.phi, psi=id.covs.psi), id.betas = list( phi=id.betas.phi, psi=id.betas.psi)
#'  , intercepts = list(phi = c(0.65,0.55,0.79,0.8), psi = 0.3 ,  r = 0.2))
#'    
#'  
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 

SimulateZ.MultiState <- function( N0 = NULL
                                , n.individuals = 100
                                , n.occasions = 5
                                , init.state
                                , state.trans
                                , plot.check = TRUE
                                , id.covs=NULL
                                , id.betas=NULL
                                , intercepts = NULL )
   {
   ##-----------------------------------
   ## ------ I. GET OBJECTS READY ------
   ##-----------------------------------
   if(!is.null(N0)){n.individuals <- sum(N0) ; n.occasions <- length(N0)}
   if(is.null(N0)){N0 <- rep(round(n.individuals/n.occasions), n.occasions)}
   z <- matrix(NA, n.individuals, n.occasions)
   f <- rep(1:length(N0), N0)
   
   ## ---- FORMAT INITIAL-STATE & STATE-TRANSITION MATRICES
   OMEGA <- array(NA, c(dim(state.trans)[1], dim(state.trans)[1], n.occasions, n.individuals))
   OMEGA[ ] <- state.trans
   
   ##-----------------------------------------------
   ## ------ II. ADD A COVARIATE ON SURVIVAL   -----
   ##-----------------------------------------------

   if(is.null(id.covs)!=TRUE){
      # names of the parameters in the list of covs, betas and intercept should similar
      if(sum(names(id.betas) %in% names(id.covs)) != length(names(id.betas)) ){
         stop("Names of the id.betas and id.covs list are not similar!")}
      
      if(sum(names(id.betas) %in% names(intercepts)) != length(names(id.betas))){
         stop("Names of the id.betas and intercepts list are not similar!")}
      # all intercept shoulbe in the state.matrix
      for(i in 1:length(intercepts)){if(sum(grepl(names(intercepts)[i], state.trans))<1){
         stop(paste("Interecept", names(intercepts)[i], " is not in present in state.trans matrix", sep="" ))}}
      
   ## ----  FOR EACH PARAMETER in the matrix 
      for(p in 1: length(intercepts)){ 
         
         # add intercepts to every occasions if the user did not allow them to vary in time
         if(length(intercepts[[p]]) < 2){ 
            intercepts[[p]] <- rep(intercepts[[p]], n.occasions)
         }
                   
         # Check if the intercept p has a covariate effect  
         if(names(intercepts)[p] %in% names(id.betas)){    
            
            temp <- formula(paste("~", paste(names(id.covs[[names(intercepts)[p]]]), collapse="+"), "-1", sep=""))
            Xmat <- model.matrix(temp, id.covs[[names(intercepts)[p]]])
            ind.fixed.effects <- Xmat %*% id.betas[[names(intercepts)[p]]]
            
         # Fill up the transition matrix for all individuals based on their covariate values 
       for(i in 1:n.individuals){
          for(t in 1:n.occasions){
            OMEGA[,,t,i] <- gsub(names(intercepts)[p]
                               , inv.logit(logit(intercepts[[p]][t]) + ind.fixed.effects[i,]) 
                               , OMEGA[,,t,i] ) 
            #                 
         }#t
       }#i
            
          # If there are no covariate effect on the parameter just replace the intercept in the transition matrix
         }else{
            for(t in 1:n.occasions){
             OMEGA[,,t,]  <- gsub(names(intercepts)[p]
                                 , intercepts[[p]][t]
                                 , OMEGA[,,t,] ) 
            }#t
          }
       }#p    
   }
    
   OMEGA <- apply(OMEGA, c(1:4), function(x) eval(parse(text=x))  )
   ##-----------------------------------------------
   ## ------ III. FORMAT INITS MATRICES ------------
   ##-----------------------------------------------
   ## ---- FORMAT INITS MATRICES
   x <- ifelse(is.vector(init.state),length(init.state), dim(init.state)[1])
   OMEGA.INIT <- array(NA,c(x, n.occasions, n.individuals))
   OMEGA.INIT[ ] <- init.state
   
   ## ---- FILL IN THE z MATRIX
   for (i in 1:n.individuals)
      {
      z[i,f[i]] <- which(rmultinom(1, 1, OMEGA.INIT[ ,f[i],i])==1)
      for (t in f[i]:(n.occasions-1))
         {
         z[i,t+1] <- which(rmultinom(1, 1, OMEGA[z[i,t], ,t,i])==1)
         }
      }
   
   
   ##-----------------------------------------------
   ## ------ IV. PLOT CHECK ------------------------
   ##-----------------------------------------------
   if(plot.check){
      if(is.null(id.covs)!=TRUE){par(mfrow=c(1,2), mar=c(5,5,5,5))}
         GetPopCompo(z, cex=0.5,lwd=0.7)
      
         # PLOT the covariate effects
         if(is.null(id.covs)!=TRUE){
            for(p in 1:length(id.covs)){
               for(ncov in 1:ncol(id.covs[[p]])){
            myX <- seq(min(id.covs[[p]][,ncov]), max(id.covs[[p]][,ncov]), length.out = 1000)
            myY <- inv.logit(logit(intercepts[[p]][1]) + id.betas[[p]][ncov] * myX)
            plot(myX, myY, xlab = names(id.covs[[p]][ncov]), ylab = names(intercepts)[[p]],
                 ylim = c(0,max(myY)), col = rev(heat.colors(length(myX))))
               }
            }
         }
      
      # PLOT intercept of each parameter for every occasions
      if(length(intercepts[[1]])>1){
         plot(1, type="n", xlab="Occasions", ylab="intercept", xlim=c(0, n.occasions+1), ylim=c(0, 1))
         col=rainbow(length(intercepts))
         for(r in 1 : length(intercepts)){
            points(intercepts[[r]]~ c(1:length(intercepts[[r]])), type="o",pch=16, col=col[r])
         }
         legend("topleft", legend=paste(names(intercepts)), lty=1, cex=0.6, col=col)
      }
   
   }
   ## ----- RETURN z MATRIX
   
   return(z)
   
}

