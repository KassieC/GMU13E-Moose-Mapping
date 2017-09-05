##In this file, we build on the results from the previous run. By examining the confusion matrix and RF model output, we identified several several groups that were at too low N to effectively code, and several related LCs that could be effectively consolidated into a common LC. Previously, in LCrecode we dropped a few catagories, and lumped a few rare catagories. Now we created a new list, LC_II_recode, that we use to train a new model, and produce our final LC list. As we went on, more things were lumped into their Vierik Level 2, 3 equivellents. 
wd="" #The directory that the results from the previous classification run (so we don't have to re-run all the rescaling/cropping).
library(raster)
library(rgdal)
library(spatstat)
library(caret)
library(randomForest)
library(e1071)
library(snow)
library(plyr)

setwd(wd)
LandSat=stack("LandSatFall_RS.tif")
solar=raster("AreaSol_RS.tif")
DEM=raster("DEM_RS.tif")
slope=raster("SLP_RS.tif")
IFSAR=raster("IFSAR_RS.tif")
LC=raster("LC_RS.tif")
SummerSat=stack("LandSatSum_RS.tif")
TRI=raster("TRI.tif")
envcov=stack(LandSat,solar,DEM,slope,IFSAR,LC,SummerSat,TRI)

names(envcov)=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "SOLAR" , "DEM", "SLP", "IFSAR", "LC", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "TRI")

LCCODEX=read.delim("recodelist.txt",header = TRUE,sep = "\t")

studyextent=readOGR(dsn = "cloudintersect_AKALB.shp", layer = "cloudintersect_AKALB")

Sys.time()
mark=Sys.time()
pointn=4000
randompts=sampleStratified(envcov$LC,size=pointn,sp=TRUE,exp=40)
dummy=as.data.frame(rep(1,nrow(randompts)))
randompts_data=SpatialPointsDataFrame(randompts,data=dummy)
dummycovs=extract(envcov,randompts)
randompts_data@data = as.data.frame(dummycovs)
mark2=Sys.time()
sampletime=mark2-mark


model_data=randompts_data@data
model_data$LC=join(randompts_data@data,LCCODEX,by="LC",type="left")$LC_II_recode ####NOTE this is the big different line here! LC_II vs. LCrecode.###
model_data=model_data[order(model_data$LC),] #This does not always work due to memory limitations. 


Sys.time()
mark=Sys.time()
modelRF = train(as.factor(LC) ~ B2 + B3 + B4 + B5 + B6 + B7 + B10 + B11 + SOLAR + DEM + SLP + IFSAR + S2 + S3 + S4 + S5 + S6 + S7 + S10 + S11, method = "rf", data=subset(model_data[complete.cases(model_data),], (LC != 16)), tuneLength = 6)
mark2=Sys.time()
modelfittime=mark2-mark
modelRF


Sys.time()
mark=Sys.time()
beginCluster()
preds_rf <- clusterR(envcov, raster::predict, args = list(model = modelRF))
endCluster()
mark2=Sys.time()
predicttime=mark2-mark

writeRaster(preds_rf,"LC_OUTPUT_30M_FS_9-22b.tif")
