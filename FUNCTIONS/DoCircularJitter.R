DoCircularJitter<-function(data.sp,r){
   
   x<-1
   orig.sp<-data.sp
   if(is.na(proj4string(orig.sp)))proj4string(data.sp)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
   
   out<-lapply(1:length(data.sp),function(x){
      this.pt<-data.sp[x,]
      circle<-raster::buffer(this.pt,r)
      out1<-spsample(circle,1,type="random")
      return(out1)
   })
   
   out<-do.call(rbind,out)
   out<-as(out,"SpatialPointsDataFrame")
   out@data<-data.frame(coordinates(out))
   
   if(is.na(proj4string(orig.sp))){
      out<-data.frame(coordinates(out))
      coordinates(out)<-out
      
      
   }
   
   # plot(e.sp)
   # plot(out,add=TRUE)
   # plot(data.sp,add=TRUE,col="red")
   
   return(out)
   
   
   
}


