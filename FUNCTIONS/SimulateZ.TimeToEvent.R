#' @title Time-dependent z state and detection probabilities.
#'
#' @description
#' \code{SimulateZ.TimeToEvent} returns p0: a \code{vector} object with the individual probabilities to be detected according to the sampling duration and pop.dynamics.
#'                      and z: a \code{vector} object with the individual state ate the end of the sampling period.
#'                      
#' @param N A \code{Numeric} object containig the number of individuals in the population at t0.
#'	@param tmax A \code{Numeric} object containing the length of the sampling period (as proportion)
#' @param annual.s A \code{Numeric} object containing the annual survival rate.
#' @param annual.r A \code{Numeric} object containing the annual reproduction rate.
#' @param annual.b A \code{Numeric} object containing the mean litter size.
#' @param annual.p A \code{Numeric} object containing the annual detection probability litter size.
#' @param plot.check A \code{Logical} object indicating whether to plot the considered scenario or not.
#'
#' @return A \code{list} object with individual p0 and z. 

SimulateZ.TimeToEvent <- function( N 
                                   , tmax  
                                   , annual.s
                                   , annual.r
                                   , annual.b
                                   , annual.p
                                   , peak.s = NULL
                                   , rate.s = NULL
                                   , peak.r = NULL
                                   , rate.r = NULL
                                   , fixed.b = FALSE
                                   , plot.check = TRUE
                                   , myCols = NULL
                                   , myLabels = NULL
                                   , y.lab = NULL)
{
   if(length(tmax) <= 1){tmax <- c(0,tmax)}
   
   ## ==== 1.ADULT SURVIVAL ====
   if(!is.null(peak.s)){
      surv <- rbinom(n = 1, size = N, prob = annual.s)
      TaD1 <- rep(2, surv)
      TaD2 <- rnorm(n = N-surv, mean = peak.s, sd = rate.s/4)
      TaD <- c(TaD1, TaD2)}else{ 
         lambda.s <- -log(annual.s)                           # Convert annual survival rate to hazard rate
         TaD <- rexp(N, lambda.s)                             # Sample individual Time at Death
         }#else                                                    
   TaB <- rep(0,length(TaD))
   
   ## ==== 2.REPRODUCTION ====
   if(!is.null(peak.r)){
      repro <- rbinom(n = 1, size = N, prob = annual.r)
      TaR1 <- rep(2,N-repro)
      TaR2 <- rnorm(n = repro, mean = peak.r, sd = rate.r/4)
      TaR <- c(TaR1,TaR2)}else{
         lambda.r <- -log(1-annual.r)                          # Convert annual survival rate to hazard rate
         TaR <- rexp(N, lambda.r)                              # Sample individual Time at Death
         }#else                        
   
   r <- ifelse(TaR <= TaD & TaR <= tmax[length(tmax)], 1, 0)       
   B <- r * rpois(length(TaR), annual.b)                        # Number of newborns produced
   if(fixed.b == TRUE){B <- r * rep(annual.b,length(TaR))} 
   TaB2 <- rep(TaR, B)                                         # Individual Time at Birth
   TaB <- c(TaB, TaB2)
   
   ## ==== 3.NEWBORN SURVIVAL ====
   if(!is.null(peak.s)){
      proba <- annual.s + pnorm(TaB, mean = peak.s, sd = rate.s/4) 
      proba <- ifelse(proba >= 1, 1, proba)
      surv2 <- rbinom(n = length(TaB2), size = 1, prob = proba)
      TaD2 <- rnorm(n = length(TaB2), mean = peak.s-TaB2, sd = rate.s/4)
      TaDD <- TaD2 + surv2}else{
         TaDD <- rexp(length(TaB2), lambda.s)           # Sample individual time alive
         }#else
   
   TaDD <- TaB2 + TaDD                                  # Calculate newborns' Time at Death
   TaD <- c(TaD,TaDD)

   ## ==== 4.DETECTION PROBABILITY ====
   lambda.p <- -log(1-annual.p)                         # Convert annual detection to hazard rate 
   Pop <- NULL
   Pop[1] <- N
   z <- p <- matrix(NA, length(TaD), length(tmax)-1)
   for(t in 2:length(tmax)){
      t0 <- apply(cbind(TaB, tmax[t-1]), 1, max)
      t1 <- apply(cbind(TaD, tmax[t]), 1, min)
      TA <- t1 - t0
      p[ ,t-1] <- 1-exp(-lambda.p*TA)
      
      z[ ,t-1] <- ifelse(TaD > tmax[t] & TaB <= tmax[t], 1, 0) # If T.a.D superior to sampling duration: z=1 ; z=0 otherwise
      Pop[t] <- sum(z[ ,t-1])
      }#t
   p[p<0] <- 0
   
   ## ==== 5.POPULATION SIZE ===
   Pop.t0 <- N
   Pop.tmax <- Pop[tmax[length(tmax)]]
   Pop.Tot <- dim(z)[1]
   
   #---------------------------------------------------------------------------------------------------------------------
   ## ==== 5.PLOT ====
   if(plot.check){
      ## ==== 1.Time to Event functions ====
      time <- seq(0,1,0.01)
      if(!is.null(peak.s)){St <- 1-(1-annual.s)*pnorm(q = time, mean = peak.s, sd = rate.s/4)}else{St <- exp(-lambda.s*time)}
      if(!is.null(peak.r)){Rt <- annual.r*pnorm(q = time, mean = peak.r, sd = rate.r/4)}else{Rt <- 1-exp(-lambda.r*time)}
      Pt <- 1-exp(-lambda.p*time)
      
      ## ==== 2.Proportion of individuals available/detected ====
      t <- seq(0,tmax[length(tmax)],0.01)
      p.available <- p.detected <- NULL
      
      # 2.1.Survival #####
      if(!is.null(peak.s)){
         surv <- rbinom(n = 1, size = 10000, prob = annual.s)
         T.a.D1 <- rep(2,surv)
         T.a.D2 <- rnorm(n = 10000-surv, mean = peak.s, sd = rate.s/4)
         T.a.D <- c(T.a.D1,T.a.D2) 
      }else{T.a.D <- rexp(10000,lambda.s)}
      
      # 2.2.Reproduction ####
      if(!is.null(peak.r)){
         repro <- rbinom(n = 1, size = 10000, prob = annual.r)
         T.a.R1 <- rep(2,10000-repro)
         T.a.R2 <- rnorm(n = repro, mean = peak.r, sd = rate.r/4)
         T.a.R <- c(T.a.R1,T.a.R2)
      }else{T.a.R <- rexp(10000, lambda.r)}                          # Sample individual Time at Death
      r <- ifelse(T.a.R <= T.a.D & T.a.R <= tmax[length(tmax)], 1, 0)              # Reproduction indicator (0 if after sampling)
      B <- r*rpois(10000, annual.b)
      T.a.B <- rep(T.a.R,B)  
      
      # 2.3.Newborn Survival ####
      if(!is.null(peak.s)){
         proba <- annual.s + pnorm(T.a.B, mean = peak.s, sd = rate.s/4) 
         proba <- ifelse(proba >= 1, 1, proba)
         surv2 <- rbinom(n =  length(T.a.B), size =1, prob = proba)
         T.a.D2 <- rnorm(n = length(T.a.B), mean = peak.s-T.a.B, sd = rate.s/4)
         T.a.DD <- T.a.D2+surv2
         }else{T.a.DD <- rexp(length(T.a.B), lambda.s)}              # Sample individual time alive
      T.a.DD <- T.a.B + T.a.DD                                       # Calculate individual Time at Death
      
      for(x in 1:length(t)){
         z <- ifelse(T.a.D >= t[x], 1, 0)
         z2 <- ifelse(T.a.B <= t[x] & T.a.DD >= t[x], 1, 0)   
         z <- c(z,z2)
         p.available[x] <- sum(z)/10000
         
         T.available <- ifelse(T.a.D > t[x], t[x], T.a.D)
         T.available2 <- ifelse(T.a.B > t[x], 0, ifelse(T.a.DD > t[x], t[x]-T.a.B, T.a.DD-T.a.B)) 
         T.available <-  c(T.available, T.available2) 
         p0 <- 1-exp(-lambda.p*T.available)
         p.detected[x] <- sum(p0)/10000
         }
      
      ## ==== 3.Plot ====
      if(is.null(myCols)){myCols <- c("firebrick3","darkgreen","navyblue","cadetblue2","cadetblue4")}
      if(is.null(y.lab)){y.lab <- seq(0,2,0.5)}
      par(mar=c(4,4,3,5)) 
      plot(time, St, type = "n", ylim = c(0,max(y.lab)), xlim = c(0,1), ylab="", xlab="", axes = FALSE)
      polygon(c(seq(0,tmax[length(tmax)],length.out=length(p.available)), rev(seq(0,tmax[length(tmax)],length.out=length(p.available)))), 
              c(p.available, rep(0,length(p.available))), border = F, col = myCols[4])
      polygon(c(seq(0,tmax[length(tmax)],length.out=length(p.detected)), rev(seq(0,tmax[length(tmax)],length.out=length(p.detected)))), 
              c(p.detected, rep(0,length(p.detected))), border = F, col = myCols[5])
      axis(4, at = y.lab, labels = y.lab*N, cex.axis = 1.2, las=1, tck=0.01, col = myCols[5], col.axis=myCols[5])
      mtext(myLabels[2], side = 4 ,col = myCols[5], line = 2.5)
      
      par(new=T) 
      plot(time, St,  ylim=c(0,1), type="l", col = myCols[1], lty=1, lwd=3, ylab="", xlab="",axes = FALSE)
      points(time, Rt, type="l", col = myCols[2], lty=1,lwd=3)
      points(time, Pt, type="l", col = myCols[3], lty=1, lwd=3)
      abline(v = tmax, lty=2)
      axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2), cex.axis = 1.2, las=1, tck=0.01)
      mtext (myLabels[1], side=2, cex=1.2, font=2)
      
      axis(1, at = seq(0,1,length.out=12), labels = seq(1,12,length.out=12), tck=0.01, cex.axis = 1.2, las = 1, hadj=0.5)
      mtext (myLabels[3], side = 1, cex = 1.2, font = 2)
      }
   
   #-----------------------------------------------------------------------------------------------------------------
   ## ==== 6.OUTPUT ====
   return(list( z = z
              , p0 = p
              , Pop = Pop
              , Pop.t0 = Pop.t0
              , Pop.tmax = Pop.tmax
              , Pop.Tot = Pop.Tot))
   }