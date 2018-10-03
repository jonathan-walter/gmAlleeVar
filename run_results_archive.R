## Get Allee threshold estimates for spatial units

rm(list=ls())

library(locfit)
library(sp)
library(raster)
library(rgdal)


setwd("/Volumes/GoogleDrive/My Drive/Gypsy moth Allee effects - 2018")

source("./RepByYear.R")
source("./AlleeParameters_flexmax.R")

## A. load spatial units and raw data
STSZones.ext<-readOGR("STSSpreadZones_extend.shp")

southeast<-read.csv("southeast_raw.csv")
uppermw<-read.csv("north_raw.csv")
midwest<-read.csv("midwest_raw.csv")

all.counts<-rbind(southeast, uppermw, midwest)
all.pts<-SpatialPointsDataFrame(coords=all.counts[,c(1,2)], data=all.counts)
proj4string(all.pts)<-proj4string(STSZones)

## Set some analysis parameters ##
bdys<-STSZones.ext #set which spatial unit boundaries to use
unique.name<-STSZones.ext$Id
nmin<-50 #minimum number of observations in a spatial unit to perform Allee parameter estimation
years<-1996:2016 #which years to get estimates
do.allyr<-TRUE
pbin<-1 #width of population bins
smpar<-0.5 #smoothing paramter for lowess
max.zt<-35 #maximum value of zt to consider

## Set up output: list of matrices. Matrices are nlocs by 3, each matrix is for a time period
output.filename="AlleePars_STSext_1996to2016.csv"
outlist<-list()

if(do.allyr==T){ #do all-year analysis if it's called for
  
  out.tt<-matrix(NA, length(bdys), 3)
  colnames(out.tt)<-c("thresh","b0","b1")
  
  for(ii in 1:length(bdys)){
    
    bdy<-bdys[ii,]
    bdycnt<-all.pts[!is.na(unlist(over(all.pts,bdy))),]
    
    if(nrow(bdycnt@data)<nmin){out.tt[ii,]<-rep(-1111,3);next}
    if(length(unique(ceiling(bdycnt@data$zt)))<10){out.tt[ii,]<-rep(-9999,3);next}
    
    bdy.RepByYear<-RepByYear(indata=bdycnt@data, pbin=pbin, smpar=smpar)
    out.tt[ii,]<-AlleeParameters_flexmax(bdy.RepByYear,year=NA, smpar=smpar)
  }
  
  outlist[["AllYear"]]<-out.tt
  
}

for(year in years){
  out.tt<-matrix(NA, length(bdys), 3)
  colnames(out.tt)<-c("thresh","b0","b1")
  
  yr.pts<-all.pts[all.pts@data$t==year,]
  
  for(ii in 1:length(bdys)){
    
    bdy<-bdys[ii,]
    bdycnt<-yr.pts[!is.na(unlist(over(yr.pts,bdy))),]
    
    if(nrow(bdycnt@data)<nmin){out.tt[ii,]<-rep(-1111,3);next}
    if(length(unique(ceiling(bdycnt@data$zt)))<10){out.tt[ii,]<-rep(-1111,3);next}
    
    bdy.RepByYear<-RepByYear(indata=bdycnt@data, pbin=pbin, smpar=smpar)
    out.tt[ii,]<-AlleeParameters_flexmax(bdy.RepByYear, year=NA, smpar=smpar)
  }
  
  outlist[[paste0("y",year)]]<-out.tt
}

## change output to matrix format
outmat<-NULL
for(ll in 1:length(outlist)){
  outmat<-cbind(outmat, outlist[[ll]])
}

hdr<-colnames(outmat)
hdr.tt<-rep(c("all","y96","y97","y98","y99","y00","y01","y02","y03","y04","y05","y06","y07","y08","y09","y10","y11",
              "y12","y13","y14","y15","y16"), each=3)
colnames(outmat)<-paste(hdr.tt, hdr, sep=".")
rownames(outmat)<-unique.name

