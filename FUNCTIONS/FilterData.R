#' @title Data set filtering
#'
#' @description
#' \code{FilterData} returns a \code{SpatialPointsDataFrame} object with the locations of the DNA samples and associated individual infos.
#'
#' @param myData A \code{SpatialPointsDataFrame} object containig cleaned NGS and dead recoveries data
#' @param poly A \code{SpatialPolygonsDataframe} with the focal study area 
#' @param country A \code{Character} with the country to subset the data to 
#' @param time_range A \code{Vector} object with the starting and ending years of interest
#' @param months A \code{Vector} object with the months to be kept. 
#' @param sex A \code{Character} object with the sex to be kept ("Hann", "Hunn","Ukjent")
#' @param dead.recovery A \code{logical} whether dead recoveries are contained in the myData and should be returned as separate sp objects.
#' @param setSex A \code{logical} whether information from all samples from a given individuals should be used to assign its sex. The function does not deal with individuals being assigned different sexex 
#' @return A \code{SpatialPointsDataFrame} object with the x and y locations of the different samples with associated attributes 
#' Id: the individual id
#' DNAID: the sample Id
#' Year: the sample biological year
#' Date: the actual dinfing date  
#' Species: the identified species
#' Country: the country the sample was found in
#' Birth. alleged birth year (when available; based on dead recoveries)
#' Death: biological year of death (for dead recoveries)
#' Age: alleged age (when available; based on dead recoveries)

FilterData <- function( myData                          ## Clean Data SpatialPoints DataFrame
                        , poly = NULL                     ## Polygon of the study area 
                        , country = NULL                  ## Alternative subset by country ("F","G","N","R","S")
                        , time_range = NULL               ## Years of interest (i.e. c(1910,2057))
                        , months = NULL                   ## Months to be selected
                        , sex = NULL                      ## "Hann", "Hunn","Ukjent"
                        , dead.recovery = FALSE
                        , setSex = FALSE){
   
   # SET THE TYPE OF PROJECTION FOR THE DATA
   proj4string(myData) <- CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
   # SET THE SEX OF INDIVIDUALS BASED ON ALL INFORMATION AVAILABLE
   if(setSex){
      # list all individual IDs
      ID <- unique(as.character(myData$Id))
      myData$Sex <- as.character(myData$Sex)
      
      # initialize the vector of IDs with conflicting sexes
      IdDoubleSex <- 0
      counter <- 1
      
      for(i in 1:length(ID)){
         # subset data to individual i
         tmp <- myData$Sex[myData$Id == ID[i]] 
         # create a table of the number of times individual i was assigned to each sex
         tab <- table(tmp[tmp %in% c("Hunn","Hann")])
         # If conflicting sexes (ID identified as both "Hunn" and "Hann")
         if(length(tab) == 2){
            # If ID assigned the same number of times to the 2 sexes, assign to Ukjent
            if(tab[1] == tab[2]){
               myData$Sex[myData$Id == ID[i]] <- "Ukjent"
            } else {
               # Otherwise pick the most common sex
               myData$Sex[myData$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
            }
            # In any case, print a warning
            print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))])) 
            IdDoubleSex[counter] <- ID[i]
            counter <- counter + 1
         }
         # If only one of "Hunn" or "Hann" registered
         if(length(tab) == 1){myData$Sex[myData$Id == ID[i]] <- names(tab)}
         # If anything else registered : "Ukjent"
         if(length(tab) == 0){myData$Sex[myData$Id == ID[i]] <- "Ukjent"}
      }#i
   }
   
   # Subset to sex
   if(!is.null(sex)){myData <- myData[which(myData$Sex %in% sex), ]} 
   
   ## Remove all samples outside the polygon of interest
   if(!is.null(poly)){myData <- myData[!is.na(over(myData, as(poly,"SpatialPolygons"))), ]}
   
   ## Alternative spatial filtering by country
   if(!is.null(country)){myData <- myData[myData$Country %in% country, ]}  
   
   ## Subset to the years of interest
   if(!is.null(time_range)){myData <- myData[which(myData$Year %in% time_range), ]} 
   
   ## Subset to the months of interest
   if(!is.null(months)){myData <- myData[which(myData$Month %in% months), ]} 
   

   if(dead.recovery){
      myData.dead <- myData[!is.na(myData$Death), ]
      myData.alive <- myData[is.na(myData$Death), ]
      myData.dead@data <- droplevels(myData.dead@data)
      myData.alive@data <- droplevels(myData.alive@data)
      IdDoubleDead <- myData.dead@data$Id[duplicated(myData.dead@data$Id)]
      
      if(setSex){
         return( list( alive = myData.alive,
                       dead.recovery = myData.dead,
                       IdDoubleSex = IdDoubleSex,
                       IdDoubleDead = IdDoubleDead))  
      }else{return( list( dead.recovery = myData.dead,
                          alive = myData.alive, 
                          IdDoubleDead = IdDoubleDead))}
   }else{
      myData@data <- droplevels(myData@data)
      return(myData)
   }
}