
#' @title Get summary from rjags output
#'
#' @description
#' \code{GetSummaryrjags} Get a summary of the rjags (coda.samples) output that also include the posteriors values for each parameter  
#'
#' @param jagsouput A \code{"mcmc.list"} object from coda.samples()functions in rjags
#' @param variable.names A \code{vector} of strings with the variable names for which gelman values should be computed 


#' @return A \code{list} object summary of the model 
#' $summary: summary of posterrios 
#' $sims.list: the list of posterior values for each parameter.  
#' @example 
#' 
#'  summary <- GetSummaryRjags( jagsouput=jagsouput
#'                            , variable.names = params )




GetSummaryRjags <- function(  jagsouput = par.samples
                              , variable.names = params){
   
   
m.mcmc <- as.matrix(jagsouput,iters = F, chains = F)

## FIND POSITION OF VARIABLES    
colnamesx <- colnames(m.mcmc)
param.pos <- NA   
for(i in 1:length(variable.names)){
   pos <- grep(variable.names[i],colnamesx)
   param.pos[(length(param.pos)+1):(length(param.pos)+length(pos))] <- pos
}
param.pos <- param.pos[!is.na(param.pos)]


## ORDER THE MATRIX WITH THE POSITION OF THE VARIABLES 
m.mcmc1 <- m.mcmc[,param.pos]

## GET SUMMARY STATISTICS 
## IN A DATAFRAME
mean <- t(t(apply(m.mcmc1,2,function(x) mean(x))))
median <- t(t(apply(m.mcmc1, 2, function(x) median(x) )))   
lch <- t(t(apply(m.mcmc1,2,function(x) quantile(x,prob=c(0.025)))))   
lcl <- t(t(apply(m.mcmc1,2,function(x) quantile(x,prob=c(0.975)))))   
df <- data.frame(cbind(mean, median, lch, lcl))  
colnames(df)  <- c("mean","median","CI2.5","CI97.5")  
df <- round(df,digit=3)


## AS A LIST FOR EACH VARIABLES 
mean <- data.frame(t(data.frame(df[,1])))
colnames(mean) <- row.names(df)
row.names(mean) <- "mean"
median <- data.frame(t(data.frame(df[,2])))
colnames(median) <- row.names(df)
row.names(median) <- "median"
CI2.5 <- data.frame(t(data.frame(df[,3])))
colnames(CI2.5) <- row.names(df)
row.names(CI2.5) <- "CI2.5"
CI97.5 <- data.frame(t(data.frame(df[,4])))
colnames(CI97.5) <- row.names(df)
row.names(CI97.5) <- "CI97.5"


## return 

return(list( summary = df
            , mean = mean
            , median = median
            , CI2.5 = CI2.5
            , CI97.5 = CI97.5))

}



