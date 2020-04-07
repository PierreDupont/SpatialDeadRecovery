#' @title Age matrix generation script for jags model with age-specific estimates
#'
#' @description
#' \code{MakeAgeMat} returns a \code{matrix} object with the age at each capture occasion 
#' for all individuals with known age at death 
#'
#' @param myData A \code{SpatialPointsDataFrame} object containig cleand up data
#' @return A \code{matrix} object with the individual age.

MakeAgeMat <- function( myData
                      , Id = NULL)                     
   {
   ## General attributes
   myData <- as.data.frame(myData)
   if(is.null(Id)) {Id <- unique(myData$Id)}
   imax <- length(Id)
   Year.min <- min(myData$Year)-1
   Year.max <- max(myData$Year) 
   tmax <- length(Year.min:Year.max)

   ## Create the age matrix
   Age.mx <- matrix(NA,imax,tmax)
   
   ## Fill in the age matrix
   for(i in 1:imax)
      {
      myData.temp <- myData[myData$Id == Id[i], ]
      D <- unique(na.omit(myData.temp$Death)) - Year.min+1
      B <- unique(na.omit(myData.temp$Birth)) - Year.min+1
      if(length(D)==1 & length(B)==0)
         {
         Age.mx[i,D] <- 99
         next
         }
      if(length(D)==1 & length(B)==1)
         {
         a1 <- 0
         if(D>=tmax){D <- tmax}
         if(B<=0){a1 <- 1-B
                  B <- 1}
         Age.mx[i,B:D] <- a1:(a1+(D-B))
         }      
      }
   return(Age.mx)
}

