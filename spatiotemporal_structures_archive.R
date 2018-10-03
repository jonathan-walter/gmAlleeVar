## Plot Allee threshold time series

rm(list=ls())

library(RColorBrewer)

setwd("/Volumes/GoogleDrive/My Drive/Gypsy moth Allee effects - 2018")

dat<-read.csv("AlleePars_STSext_1996to2016.csv", na.strings = c("-1111","-5555"))
years<-1996:2016

thresh.series<-dat[,seq(5,ncol(dat),by=3)]
#thresh.series[thresh.series==-999 | thresh.series==-555]<-300

colors<-brewer.pal(12, "Set3")
colors[2]<-rgb(t(col2rgb("cadetblue")),maxColorValue = 250)

thresh.series.max<-as.matrix(abs(thresh.series))
thresh.series.na<-as.matrix(thresh.series)
thresh.series.na[thresh.series.na < 0]<-NA


## Look at spatial synchrony in AT ################

## Remove Region 12 from SS analysis!!!!!!!!!!!!!!!!!!!!!!!!!!!

cor.thresh.na<-cor(t(thresh.series.na[-12,]), method="spearman", use="pairwise.complete.obs")
print(cor.thresh.na)
summary(c(cor.thresh.na[lower.tri(cor.thresh.na)]))


cor.thresh.max<-cor(t(thresh.series.max[-12,]), method="spearman", use="pairwise.complete.obs")
print(cor.thresh.max)
summary(c(cor.thresh.max[lower.tri(cor.thresh.max)]))
# image(cor.thresh)
# print(mean(cor.thresh[lower.tri(cor.thresh)]))

dd<-c(1:10, 1:9, 1:8, 1:7, 1:6, 1:5, 1:4, 1:3, 1:2, 1)

spatSynch.na.long<-data.frame(dist=dd, corr=c(cor.thresh.na[lower.tri(cor.thresh.na)]))
spatSynch.max.long<-data.frame(dist=dd, corr=c(cor.thresh.max[lower.tri(cor.thresh.max)]))

spatSynch.max.long2<-spatSynch.max.long[spatSynch.max.long$dist <=6,]

tiff("Fig2_timeseriesSpatSynch_20180928.tif", units="in", width=6.5, height=8, res=300)
par(mfrow=c(4,1), mar=c(3.5,3.5,2.1,1), mgp=c(2.1,0.7,0), tcl=-0.4)

plot(NA,NA, xlim=c(min(years),max(years)), ylim=c(0,4), xlab="Year", ylab="log10(Allee threshold+1)",
     main="Wisconsin & Minnesota", cex.axis=1, cex.lab=1.1, cex.main=1.1, lwd=2)
for(ii in 1:4){
  lines(years, log10(thresh.series.max[ii,]+1), col=colors[ii], type="b", lwd=2)
  if(sum(is.na(thresh.series.na[ii,]))!=0){
    na.yrs<-years[which(is.na(thresh.series.na[ii,]))]
    points(na.yrs, log10(thresh.series.max[ii,][which(is.na(thresh.series.na[ii,]))]+1), col=colors[ii], pch="x", cex=2)
  }
}
legend("top", legend=paste0("R",1:4), lwd=1.5, col=colors[1:4], ncol=4, cex=0.9)
text(1995.5,3.9,"a)")

plot(NA,NA, xlim=c(min(years),max(years)), ylim=c(0,4), xlab="Year", ylab="log10(Allee threshold+1)",
     main="Ohio, Indiana, & Illinois", cex.axis=1, cex.lab=1.1, cex.main=1.1, lwd=2)
for(ii in 5:8){
  lines(years, log10(thresh.series.max[ii,]+1), col=colors[ii], type="b", lwd=2)
  if(sum(is.na(thresh.series.na[ii,]))!=0){
    na.yrs<-years[which(is.na(thresh.series.na[ii,]))]
    points(na.yrs, log10(thresh.series.max[ii,][which(is.na(thresh.series.na[ii,]))]+1), col=colors[ii], pch="x", cex=2)
  }
}
legend("top", legend=paste0("R",5:8), lwd=1.5, col=colors[5:8], ncol=4, cex=0.9)
text(1995.5,3.9,"b)")

plot(NA,NA, xlim=c(min(years),max(years)), ylim=c(0,4), xlab="Year", ylab="log10(Allee threshold+1)",
     main="Virginia, West Virginia, & North Carolina", cex.axis=1, cex.lab=1.1, cex.main=1.1, lwd=2)
for(ii in 9:11){
  lines(years, log10(thresh.series.max[ii,]+1), col=colors[ii], type="b", lwd=2)
  if(sum(is.na(thresh.series.na[ii,]))!=0){
    na.yrs<-years[which(is.na(thresh.series.na[ii,]))]
    points(na.yrs, log10(thresh.series.max[ii,][which(is.na(thresh.series.na[ii,]))]+1), col=colors[ii], pch="x", cex=2)
  }
}
legend("top", legend=paste0("R",9:11), lwd=1.5, col=colors[9:11], ncol=3, cex=0.9)
text(1995.5,3.9,"c)")

boxplot(corr~dist, data=spatSynch.max.long2, main="Spatial synchrony", ylab="Correlation", xlab="Distance",
        cex.axis=1, cex.lab=1.1, cex.main=1.1)
text(0.35,0.56,"d)")
dev.off()

##Look at temporal structure in AT ################

lagCorr<-function(zt, nlags=5, type="spearman"){
  
  out<-rep(NA, nlags)
  names(out)<-paste0("lag",1:nlags)
  
  for(lag in 1:nlags){
    
    lagzt<-zt[-c(1:lag)]
    tmpzt<-zt[1:length(lagzt)]
    out[lag]<-cor(tmpzt,lagzt,method=type, use="pairwise.complete.obs")
    
  }
  return(out)
}

lagCorr.na<-t(apply(thresh.series.na, 1, lagCorr, nlags=10))
lagCorr.max<-t(apply(thresh.series.max, 1, lagCorr, nlags=10))

lagCorr.na.long<-data.frame(region=rep(1:12,10),lag=rep(1:10,each=12),corr=c(lagCorr.na))
lagCorr.max.long<-data.frame(region=rep(1:12,10),lag=rep(1:10,each=12),corr=c(lagCorr.max))

