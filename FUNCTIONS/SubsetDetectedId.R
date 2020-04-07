

SubsetDetectedId <- function(data, detections, replace.NA = FALSE)
   {
   ## Identify detected individuals for each data source
   if(is.list(detections)){
      detected <- list()
      for(i in 1:length(detections)){
         max.detections <- max(detections[[i]])
         detected[[i]] <- apply(detections[[i]], 1, sum)>0
         if(!any(detections[[i]]==0)){detected[[i]] <- apply(detections[[i]], 1, function(x){any(x != max.detections)})}
         }#i
      detected.all <- unique(unlist(lapply(detected, which)))
   }else{
      max.detections <- max(detections)
      detected.all <- apply(detections, 1, sum)>0
      if(!any(detections==0)){detected.all <- apply(detections, 1, function(x){any(x != max.detections)})}}

   ## Subset data to detected individuals
   if(is.list(data) & !is.data.frame(data)){
      if(!replace.NA){
         for(i in 1:length(data)){
            if(is.vector(data[[i]])){data[[i]] <- data[[i]][detected.all]}
            if(length(dim(data[[i]]))==2){data[[i]] <- data[[i]][detected.all, ]}
            if(length(dim(data[[i]]))==3){data[[i]] <- data[[i]][detected.all, , ]}
            if(length(dim(data[[i]]))==4){data[[i]] <- data[[i]][detected.all, , , ]}
            if(length(dim(data[[i]]))==5){data[[i]] <- data[[i]][detected.all, , , , ]}
            if(length(dim(data[[i]]))>=6){print("An array with more than 5 dimensions??? Seriously????")}}
      }else{
         for(i in 1:length(data)){
            if(is.vector(data[[i]])){data[[i]][detected.all] <- NA}
            if(length(dim(data[[i]]))==2){data[[i]][detected.all, ] <- NA}
            if(length(dim(data[[i]]))==3){data[[i]][detected.all, , ] <- NA}
            if(length(dim(data[[i]]))==4){data[[i]][detected.all, , , ] <- NA}
            if(length(dim(data[[i]]))==5){data[[i]][detected.all, , , , ] <- NA}
            if(length(dim(data[[i]]))>=6){print("An array with more than 5 dimensions??? Seriously????")}}}
   }else{
      if(!replace.NA){
         if(is.vector(data)){data <- data[detected.all]}
         if(length(dim(data))==2){data <- data[detected.all, ]}
         if(length(dim(data))==3){data <- data[detected.all, , ]}
         if(length(dim(data))==4){data <- data[detected.all, , , ]}
         if(length(dim(data))==5){data <- data[detected.all, , , , ]}
         if(length(dim(data))>=6){print("An array with more than 5 dimensions??? Seriously????")}
      }else{
         if(is.vector(data)){data[[i]][detected.all] <- NA}
         if(length(dim(data))==2){data[detected.all, ] <- NA}
         if(length(dim(data))==3){data[detected.all, , ] <- NA}
         if(length(dim(data))==4){data[detected.all, , , ] <- NA}
         if(length(dim(data))==5){data[detected.all, , , , ] <- NA}
         if(length(dim(data))>=6){print("An array with more than 5 dimensions??? Seriously????")}}}
      
   return(data)
   }
