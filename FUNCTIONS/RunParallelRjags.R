#' @title Parallel MCMC sampling using rjags.
#'
#' @description
#' \code{RunParallelRjags} returns a \code{list} object with the results of MCMC sampling using rjags
#' either as an MCMC list (raw rjags output)
#' or as a list with various features ( summary, DIC etc...).
#' Not sure it works on Mac/linux... (if that is the case, use JagsUI parrallel function)
#'
#' @param jags.input A \code{list} object containig the data needed to run the JAGS model
#'	@param inits A \code{list} object containing the initial values for each core/chain
#' @param params A \code{vector} object containing the names of the parameters to be returned
#' @param model.file A \code{Characters} object with the name of the JAGS model 
#' @param n.cores A \code{Numeric} with the number of cores to be run in parallel (usually equal to the number of chains)
#' @param n.adapt A \code{Numeric} with the number of itertions to be run as adaptive phase for each core/chain
#' @param n.iter A \code{Numeric} with the number of iterations to be run for sampling for each core/chain
#' @param n.thin A \code{Numeric} the number of samples to be ignored between two stored samples
#'
#' @param friendly.output A \code{logical}  Whether results should be presented as raw MCMC lists (rjags) or user-friedly output (jagsUI)
#' 
#' @return A \code{list} object with the MCMC samples for each parameter of the model specified in "params" 
#' and the time needed to run each phase of the model fitting 

RunParallelRjags <- function( jags.input
                            , inits
                            , params
                            , model.file
                            , n.cores
                            , n.iter
                            , n.adapt
                            , n.thin
                            , friendly.output = TRUE
                            , DIC = FALSE
                            , params.omit = NULL)
   {
   inits1 <- inits
   inits <- list()
   for(i in 1:n.cores){
      inits[[i]] <- inits1()
      }
   
   cl <- makeCluster(n.cores, "SOCK")
   
   clusterExport(cl, c("jags.input", "inits", "params", "model.file", "n.iter", "n.adapt", "n.thin"), envir = environment(NULL))

   clusterEvalQ(cl, library(rjags))

   coda.samples.wrapper <- function(x){
      t0 <- proc.time()
      inits2  <- inits[[x]]
      inits2$.RNG.name <- "base::Wichmann-Hill"## random seed generator
      inits2$.RNG.seed <- x
      temp.model <- jags.model(file = model.file, data = jags.input, inits = inits2, n.chains = 1, n.adapt = n.adapt)
      t1 <- proc.time()
      mySamples <- coda.samples(model = temp.model, variable.names = params, n.iter = n.iter, thin = n.thin)
      t2 <- proc.time()
      adaptive.runtime <- t1-t0
      sampling.runtime <- t2-t1
      total.runtime <- t2-t0
      return(list( mySamples = mySamples
                 , adaptive.runtime = adaptive.runtime
                 , sampling.runtime = sampling.runtime
                 , total.runtime = total.runtime))
      }
   
   print(paste("Beginning parallel processing using", n.cores, "cores. Console output will be suppressed."))
   mySamples <- clusterApply(cl, 1:n.cores, coda.samples.wrapper)
   stopCluster(cl)
   
   adaptive.runtime <- sampling.runtime <- total.runtime <- list()
   for(i in 1:length(mySamples)){
      adaptive.runtime[[i]] <- mySamples[[i]]$adaptive.runtime
      sampling.runtime[[i]] <- mySamples[[i]]$sampling.runtime
      total.runtime[[i]] <- mySamples[[i]]$total.runtime
      mySamples[[i]] <- mySamples[[i]]$mySamples[[1]]
      }
   
   class(mySamples) <- "mcmc.list"
   
   if(!friendly.output){
      return(list( samples = mySamples
                 , adaptive.runtime = adaptive.runtime
                 , sampling.runtime = sampling.runtime
                 , total.runtime = total.runtime)) 
      }
   else{
      myResults <- ProcessCodaOutput( x = mySamples, DIC = DIC, params.omit = params.omit, verbose = TRUE)
      return(list( samples = mySamples
                 , JAGS.output = myResults
                 , adaptive.runtime = adaptive.runtime
                 , sampling.runtime = sampling.runtime
                 , total.runtime = total.runtime)) 
      }
   }