#' @title Wrapper function to simulate and create jags input data from a set of SCR parameters. The output is directly exported as an RData file.
#'
#' @description
#' \code{MakeJagsSimData} returns a list object with \code{} 
#' 
#' @param parms A list of named parameters (single set). Can also be a dataframe with a single row.
#' @param habitat.ls List object, output from the MakeHabitat function.
#' @param cov.RasterStack Raster stack with one raster for each covariate (standardized), in the same projection as the focal area polygon in habitat.ls.
#' @param coefs A numeric vector denoting the coefficients associated with the covariates (from the cov.RasterStack). For now the same coefficients and covariates are applied to AC placement and detection prob.
#' @param outpath A string denoting the path to the directory into which the output is to be saved.
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeJagsSimData.R
#' @keywords simul
#'
#' @examples
#' # Simulate and create jags input files.
#'MakeJagsSimData(
#'   parms=parm.df[1,],
#'   habitat.ls=habitat.ls,
#'   cov.RasterStack=cov.RasterStack,
#'   coefs=coefs,
#'   outpath="C:\\JagsInputFiles\\",
#'   plot.check = TRUE







MakeJagsSimData<-function(parms,habitat.ls,cov.RasterStack,coefs,outfilename="test.RData",plot.check=FALSE){
   
   
   if(plot.check){par(mfrow=c(2,3),mar=c(1,1,1,1))}
   
   detector.ls<-MakeSearchGrid(
      data=habitat.ls$habitat.poly,
      resolution=parms$detector.cell.size,
      div=1,
      center=TRUE,
      plot=plot.check
   )
   
   
   detector.ls$detector.scaled.xy<-UTMToGrid(
      data.sp=detector.ls$detector.sp,
      grid.sp=habitat.ls$habitat.sp,
      plot.check = plot.check)$data.scaled.xy
   
   
   
   
   #----AC PLACEMENT (FOR NOW: same betas and covariates as for detection, just associated with habitat raster cell centers)
   i<-2
   temp<-lapply(1:length(coefs),function(i){
      this.coef.name<-names(coefs)[i]
      this.coef<-coefs[i]
      this.raster<-cov.RasterStack[[this.coef.name]]
      out<-raster::extract(this.raster,habitat.ls$habitat.sp)
   })
   covs<-do.call(data.frame,temp)
   names(covs)<-names(coefs)
   betas<-c(intercept=log(parms$lambda0),coefs)
   
   #----AC PLACEMENT
   
   AC.sp<-SimulateACs(
      N=parms$n.true,
      habitat.sp=habitat.ls$habitat.sp[habitat.ls$habitat.index,]
      ,covs=covs[habitat.ls$habitat.index,],
      betas=betas,
      plot=plot.check
      )

   #----DETECTION (SETUP SPECIFICALLY FOR UNGULATE STUDY: will need to adjust for other analyses)
   i<-2
   temp<-lapply(1:length(coefs),function(i){
      this.coef.name<-names(coefs)[i]
      this.coef<-coefs[i]
      this.raster<-cov.RasterStack[[this.coef.name]]
      out<-raster::extract(this.raster,detector.ls$detector.sp)
   })
   covs<-do.call(data.frame,temp)
   names(covs)<-names(coefs)
   betas<-c(intercept=log(parms$lambda0),coefs)
   
   #covs<-data.frame(clc=covs$clc)#--since there may be a problem with the elevation cov (not enough detectors with "different" suitability values)
   
   n.covs.lambda<-length(betas)#dim(X.lambda)[2]
   
   det.sim<-SimulateDetection(p0=parms$lambda0,
                              sigma=parms$sigma,#sigma,
                              AC.sp=AC.sp,
                              detector.sp=detector.ls$detector.sp,
                              detector.covs=covs, 
                              betas=betas,
                              n.samples=parms$n.samples.total,
                              type="Poisson",#
                              plot=plot.check)
   

   ####### DATA AUGMENTATION
   
   n.individuals<-dim(det.sim$y)[1]
   n.detectors<-dim(det.sim$y)[2]
   
   #need to make more individuals available than actual pop size, hence:
   n.aug.ind<-ceiling(parms$n.true*2)
   #.... instead of:
   #n.aug.ind<-n.individuals*0.75
   
   y<-rbind(det.sim$y,array(0,dim=c(n.aug.ind-n.individuals,n.detectors)))
   
   n.individuals<-dim(y)[1]
   
   
   #--need to keep the original names for detected individuals to be able to link them with their true ACs
   dimnames(y)[[1]]<-ifelse(dimnames(y)[[1]]=="",paste("aug",1:sum(dimnames(y)[[1]]=="")),dimnames(y)[[1]])
   
   z<-as.numeric(apply(y,1,max)>0)# KNOWN STATES
   z[z==0]<-NA# UNKNOWN STATES
   
   
   
   ####### COMPILE HABITAT mx and xy
   
   habitat.xy<-habitat.ls$habitat.scaled.xy
   habitat.xy<-habitat.xy[habitat.ls$habitat.index,]
   
   habitat.mx<-habitat.ls$habitat.mx
   habitat.mx[habitat.mx>0]<-1
   
   y.max<-dim(habitat.mx)[1]#max(habitat.xy[,"y"])
   x.max<-dim(habitat.mx)[2]#max(habitat.xy[,"x"])
   
   
   ####### INITIAL VALUES
   
   z.init<-ifelse(!is.na(z),NA,rbinom(length(which(is.na(z))),1,0.5))
   sx.init<-2
   
   inits <- function(){#---consider adding starting values for sx and sy
      list("z"=z.init)}#,"alpha0"=runif(1,-1.5,-1),"sigma"=runif(1,0.1,5))}
   
   
   
   OK<-rep(1,n.individuals)
   
   
   my.jags.input <- list(X.p=det.sim$X,n.covs.p=n.covs.lambda,OK=OK,z=z,y=y, n.detectors=n.detectors, n.individuals=n.individuals, detector.xy=detector.ls$detector.scaled.xy,habitat.mx=habitat.mx ,x.max=x.max,y.max=y.max)
   #parameters<-c("N","sigma","beta.p","alpha0")#"sx","sy","z",
   
   print(paste("Saving jags input file for parameter set nr. ",parms$set.id,sep=""))
   #fileName<-paste(outpath,"JagsInputSetID",parms$set.id ,".RData",sep="")
   #this.jags.input.set<-save(parms,my.jags.input,inits,z.init,AC.sp,detector.ls,habitat.ls,file=outfilename)
   this.jags.input.set<-save(parms,my.jags.input,inits,z.init,AC.sp,file=outfilename)
}
