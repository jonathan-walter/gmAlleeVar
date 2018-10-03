## function for computing fractions of populations replacing themselves in a given year

## Code by Patrick Tobin and Jonathan Walter, last edit 2018-01-27
## contact: jaw3es@virginia.edu

RepByYear<-function(indata,pbin=1,smpar=0.5){
  #indata is a data frame containing:
    #t = year t
    #zt = the estimated population density in year t
    #zt1 = the estimated population density in year t+1
  #pbin = the width of population bins (defaults to 1)
  #smpar = the smoothing parameter for the spline (defaults to 0.5)
  
  out<-NULL
  
  years<-unique(indata$t)
  
  for(yy in years){
    
    subdat<-indata[indata$t==yy,]
    bivar<-subdat$zt < subdat$zt1
    nbin<-ceiling(subdat$zt/pbin)
    
    #out.yy<-data.frame(year=rep(yy,nrow(subdat)))
    out.yy<-data.frame(zt=sapply(split(subdat$zt,nbin), mean))
    out.yy$p.rep<-sapply(split(bivar,nbin), mean)
    
    fit.yy<-lowess(out.yy$zt, out.yy$p.rep, f=smpar, iter=0)
    out.yy$p.rep.sm<-fit.yy$y
    out.yy$year<-rep(yy, nrow(out.yy))
    out<-rbind(out, out.yy)
  }
  return(out)
}
  