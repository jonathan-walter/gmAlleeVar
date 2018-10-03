## Function for computing Allee parameters from output of RepByYear

## Code by Patrick Tobin and Jonathan Walter, last edit 2018-06-206
## contact: jaw3es@virginia.edu

## Note that this also calculates the Allee slope (Walter et al. [2017] Population Ecology),
## but it is not used further because simulation studies showed that the slope of the 
## density-replacement rate relationship was a poor estimate of the density-growth rate
## relationship (J. Walter, unpublished data)

AlleeParameters_flexmax<-function(indata, year=NA, smpar=0.5, slope.method="lin.interp", start.ztmax = 30){
  #indata is output from RepByYear
  #year = the year to estimate for, if year = NA, use all years
  #smpar =  the LOWESS smoothing parameter (defaults to 0.5)
  
  # if(slope.method != "lin.interp"){stop("only slope method == lin.interp is implemented")}
  
  library(locfit)
  
  out<-rep(NA, 3)
  names(out)<-c("thresh","b0","b1")
  

  if(!is.na(year)){
    indata<-indata[indata$t==year,]
  }

    if(max(indata$zt)<(start.ztmax+10)){maxseq=start.ztmax}
    else{maxseq<-seq(start.ztmax, max(indata$zt), by=10)}
    
    for(max.zt in maxseq){
      subdata<-indata[indata$zt<=max.zt,]
    
      fit<-locfit(subdata$p.rep~subdata$zt, deg=1, alpha=smpar)
      der<-locfit(subdata$p.rep~subdata$zt, deg=1, alpha=smpar, der=1)
      fit.pred<-predict(fit,0:max.zt)
    
    #compute Allee threshold

    if(max(fit.pred)<0.5){out<-rep(-max.zt,3);next}
    else if(fit.pred[1]>=0.5){out[1]<-0}
    else{
      nbin<-0:max.zt
      min.rep<-which.min(fit.pred<0.5)    #identify two points on either side of 0.5
      nx<-nbin[c(min.rep-1,min.rep)]
      yx<-fit.pred[c(min.rep-1,min.rep)]
      coeffs<-lm(yx~nx)$coefficients    #linear interpolation to estimate density
      out[1]<-(0.5-coeffs[1])/coeffs[2]
    }    
    #get intercept and slope
    out[2]<-predict(fit,0)
    
    if(slope.method=="lin.interp"){out[3]<-mean(predict(der,0:which.max(fit.pred)))}
    if(slope.method=="d(0)"){out[3]<-predict(der,0)}
    break
    }
  #}
    
return(out)    

}