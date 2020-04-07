#' @title Function create NIMBLE input data when some SXY posisions are "forced" due to long-distance dispersal.
#'
#' @description
#' \code{plot.violins} creates a violin plot from a list of numerical vectors.
#' 
#' @param dat.list A list of numerical vectors. Each list item will be associated with its own violin, representing its distribution.
#' @param x A character or numerical vector indicating the labels associated with each item in dat.list.
#' @param at The location of each violin on the X axis.
#' @param add If the violins to be added to an existing plot (TRUE) or if a new plot is to be generated (FALSE)
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/ForceDisperserSXY.R
#' @keywords simul
#'
#' @examples
#' # Generate a violin plot
#' ForceDisperserSXY(sxy=sxy,y.alive=y.alive,d.threshold=10,habitat.sp=myHabitat$habitat.sp,plot.check=TRUE)
#' 

ForceDisperserSXY<-function(sxy,y.alive,d.threshold,habitat.sp,plot.check=TRUE){
   
   sxy.temp<-sxy
   y.detected<-apply(y.alive,MAR=c(1,3),max)
   sxy.temp[,1,][!y.detected]<-NA
   sxy.temp[,2,][!y.detected]<-NA
   
   
   
   
   ####---DETECT DISPERSAL EVENTS
   i<-3
   out2<-lapply(1:dim(y.alive)[1],function(i){
      x<-sxy.temp[i,,]
      ii<-5
      out1<-lapply(2:dim(x)[2],function(ii){
         out<-data.frame(ind=NA,time=NA,dist=NA)
         this.x<-x[,ii]
         if(any(!is.na(x[1,1:(ii-1)]))){
         previous.x<-x[,max(which(!is.na(x[1,1:(ii-1)])))]
         if(!any(is.na(this.x)) & all(!is.na(previous.x))){
            
            this.dist<-sqrt(sum((this.x-previous.x)^2))
            #out<-this.dist
            if(this.dist>d.threshold)out<-data.frame(ind=i,time=ii,dist=this.dist)
            return(out)
         }}})
      out1<-do.call(rbind,out1)
      return(out1)
      
   })
   
   out2<- out2[!unlist(lapply(out2,is.null))]
   out2<-na.omit(do.call(rbind,out2))

   

   if(dim(out2)[1]>0){
   sxy.dispersers<-sxy.temp
   
   sxy.dispersers[]<-NA
   
   for(i in 1:dim(out2)[1]){
      
      sxy.dispersers[out2[i,"ind"],c(1,2),out2[i,"time"][[1]]]<-sxy[out2[i,"ind"],c(1,2),out2[i,"time"][[1]]]
   }
   
   
 
   
   
   mask<-sxy.dispersers
   sxy.dispersers[is.na(sxy.dispersers)]<-0
   sxy.data<-sxy.data.scaled<-sxy.dispersers
   for( t in 1:dim(sxy)[3]){
      sxy.data.scaled[,,t] <- UTMToGrid(grid.sp = habitat.sp,
                                       data.sp = SpatialPoints(sxy.dispersers[,,t]),
                                       plot.check = F
      )$data.scaled.xy
   }
   sxy.data[is.na(mask)]<-NA
   sxy.data.scaled[is.na(mask)]<-NA
   
   sxy.init<-sxy.init.scaled<-sxy
   for( t in 1:dim(sxy)[3]){
      sxy.init.scaled[,,t] <- UTMToGrid(grid.sp = habitat.sp,
                            data.sp = SpatialPoints(sxy[,,t]),
                            plot.check = F
      )$data.scaled.xy
   }
   sxy.init[!is.na(mask)]<-NA
   sxy.init.scaled[!is.na(mask)]<-NA
   
   
   if(plot.check){
      
      all.sxy<-do.call(rbind,lapply(1:dim(sxy)[3],function(x)sxy[,,x]))
      plot(all.sxy[,1],all.sxy[,2], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
      col <- rainbow(dim(sxy)[1])
      #points(sxy[,1,1],sxy[,2,1], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
      points(sxy.temp[,1,1],sxy.temp[,2,1], pch=21, bg=col,cex=1.3,col="black")
      for(t in 2:dim(sxy)[3]){
         #points(sxy[,1,t],sxy[,2,t], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
         #points(sxy.temp[,1,t],sxy.temp[,2,t], pch=21, bg=col,cex=1.3,col="black")
         arrows( x0 = sxy[,1,t-1]
                 , x1 = sxy[,1,t]
                 , y0 = sxy[,2,t-1]
                 , y1 = sxy[,2,t], col = adjustcolor(col,alpha=0.3), length = 0.01,lty=3)
      }#t
      
    
         for(i in 1:dim(sxy.temp)[1]){
            x<-sxy.temp[i,,]
             x2<-na.omit(t(as.matrix(x)))
         if(dim(x2)[1]>1){
            
            for(t in 2:dim(x2)[1]){
              arrows( x0 = x2[t-1,1]
                       , x1 = x2[t,1]
                       , y0 = x2[t-1,2]
                       , y1 = x2[t,2], col = col[i], length = 0.08,lwd=2)
            }#t
            
         }
         
      }
      
      
      apply(out2,1,function(x){
         points(sxy[x["ind"][[1]],1,x["time"][[1]]],sxy[x["ind"][[1]],2,x["time"][[1]]],cex=2,pch=6,lwd=2,col="black")
      })
      
      # apply(sxy.dispersers,3,function(x){
      #    points(x[,1],x[,2],pch=4,col="red",lwd=2,cex=1.5)
      #    
      # })
      
   }
   }
   if(dim(out2)[1]==0){
      
      sxy.init.scaled<-sxy.init<-sxy.data<-sxy.data.scaled<-sxy
      sxy.data[]<-sxy.data.scaled<-NA

      for( t in 1:dim(sxy)[3]){
         sxy.init.scaled[,,t] <- UTMToGrid(grid.sp = habitat.sp,
                                           data.sp = SpatialPoints(sxy.init[,,t]),
                                           plot.check = F
         )$data.scaled.xy
      }
      
      if(plot.check){
         
         all.sxy<-do.call(rbind,lapply(1:dim(sxy)[3],function(x)sxy[,,x]))
         plot(all.sxy[,1],all.sxy[,2], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
         col <- rainbow(dim(sxy)[1])
         #points(sxy[,1,1],sxy[,2,1], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
         points(sxy.temp[,1,1],sxy.temp[,2,1], pch=21, bg=col,cex=1.3,col="black")
         for(t in 2:dim(sxy)[3]){
            #points(sxy[,1,t],sxy[,2,t], pch=21, bg=adjustcolor(col,alpha=0.3),cex=1,col="white")
            #points(sxy.temp[,1,t],sxy.temp[,2,t], pch=21, bg=col,cex=1.3,col="black")
            arrows( x0 = sxy[,1,t-1]
                    , x1 = sxy[,1,t]
                    , y0 = sxy[,2,t-1]
                    , y1 = sxy[,2,t], col = adjustcolor(col,alpha=0.3), length = 0.01,lty=3)
         }#t
         
         
         for(i in 1:dim(sxy.temp)[1]){
            x<-sxy.temp[i,,]
            x2<-na.omit(t(as.matrix(x)))
            if(dim(x2)[1]>1){
               
               for(t in 2:dim(x2)[1]){
                  arrows( x0 = x2[t-1,1]
                          , x1 = x2[t,1]
                          , y0 = x2[t-1,2]
                          , y1 = x2[t,2], col = col[i], length = 0.08,lwd=2)
               }#t
               
            }
            
         }
         

         
      }
      
   }
   out<-list(sxy.init=sxy.init,sxy.data=sxy.data,sxy.init.scaled=sxy.init.scaled,sxy.data.scaled=sxy.data.scaled)
   return(out)
   
}

