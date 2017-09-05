####This code was origionally written for the SU-WA project, to take the AK NHP/ACCS data (western AK LC mosiac) and convert the study area to a consistently applied set of classifications. Our approach was to take imagery and ancillary data (DEM, TRI, so on) and use the AK NHP LC data as 'points of truth' (with a small t) and create a random forest scheme to predict those points from the imgagery/predictive layers. We then take the enviromental covariates and create a new predicted Land Cover raster. 

sourcedir= #The source directory for the imagery
wd=  #Where we will eventually output product to.
setwd(sourcedir)
library(raster)
library(rgdal)
library(spatstat)
library(caret)
library(randomForest)
library(e1071)
library(snow)
library(plyr)

LandSat=stack("LS80690162014268_AKALB.tif") #Landsat imagery

AreaSol=raster("AreaSol.tif") #Solar radiation for the scene

DEM=raster("DEMClp.tif") #NED DEM for the area
converttometers<- function(x) {x*0.3048} #Write a brief function to convert ft to meters
temp=calc(DEM,converttometers) #Make metric
DEM=temp  #Overwrite
rm(temp) #Remove the intermediate file to free up memory (memory management is going to be important here).

# terrain(DEM,opt="TRI",filename="TRI.tif") #The origional line for generating the TRI raster
TRI=raster("TRI.tif")

IFSAR=raster("ifsarClp_AKALB.tif") #Using GINA IFSAR data, we took DSM-DTM, which gave us how much the DSM varied over the terrain. We consider this an index of canopy height.

LC=raster("NHPVegClp_AKALB.tif") #The AK Natural Heritage Program's LC converted to AK Albers

SLP=raster("SLPClp.tif") #Slope

SummerSat=stack("LC80690162016242_AKALB.tif") #Landsat imagery for a 2nd day.

#Now we create a common bounding box for all the image sets so we have a common raster extent. 
xmin <- max(bbox(LandSat)[1,1], bbox(AreaSol)[1,1], bbox(DEM)[1,1],bbox(IFSAR)[1,1],bbox(LC)[1,1],bbox(SLP)[1,1],bbox(SummerSat)[1,1])
xmax <- min(bbox(LandSat)[1,2], bbox(AreaSol)[1,2], bbox(DEM)[1,2],bbox(IFSAR)[1,2],bbox(LC)[1,2],bbox(SLP)[1,2],bbox(SummerSat)[1,2])
ymin <- max(bbox(LandSat)[2,1], bbox(AreaSol)[2,1], bbox(DEM)[2,1],bbox(IFSAR)[2,1],bbox(LC)[2,1],bbox(SLP)[2,1],bbox(SummerSat)[2,1])
ymax <- min(bbox(LandSat)[2,2], bbox(AreaSol)[2,2], bbox(DEM)[2,2],bbox(IFSAR)[2,2],bbox(LC)[2,2],bbox(SLP)[2,2],bbox(SummerSat)[2,2])

#Write new extent and crop all the rasters. This code is slow.
newextent=c(xmin,xmax,ymin,ymax)

LandSat=crop(LandSat,newextent)
SummerSat=crop(SummerSat,newextent)
AreaSol=crop(AreaSol,newextent)
DEM=crop(DEM,newextent)
TRI=crop(TRI,newextent)
IFSAR=crop(IFSAR,newextent)
LC=crop(LC,newextent)
SLP=crop(SLP,newextent)

#Resample the cropped rasters to a common grid - using the AK NHP LC as the common grid. This code is glacially slow.
LandSat=resample(LandSat,LC)
AreaSol=resample(AreaSol,LC)
DEM=resample(DEM,LC)
IFSAR=resample(IFSAR,LC) 
#LC=resample(LC,LC)## Obviously, this line makes no sense, but is included for symmetry reasons.
SLP=resample(SLP,LC)
SummerSat=resample(SummerSat,LC)
TRI=resample(TRI,LC)

#Write the outputs to disk so you don't have to re-run the code. Note that _RS files are all 'resampled.'
setwd(wd)
writeRaster(LandSat,"LandSatFall_RS.tif")
writeRaster(SummerSat,"LandSatSum_RS.tif")
writeRaster(AreaSol,"AreaSol_RS.tif")
writeRaster(DEM,"DEM_RS.tif")
writeRaster(SLP,"SLP_RS.tif")
writeRaster(IFSAR,"IFSAR_RS.tif")
writeRaster(LC,"LC_RS.tif")
writeRaster(TRI,"TRI.tif")

#######Checkpoint#########

#Load the RS rasters from disk and stack them into a brick.
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

names(envcov)=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "SOLAR" , "DEM", "SLP", "IFSAR", "LC", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "TRI") #B1-11 is band 1-11 for the first imgagery set. S1-11 is Band 1-11 for the second. 

LCCODEX=read.delim("recodelist.txt",header = TRUE,sep = "\t") #This is a list of integer/catigorical data and the text string land cover type they correspond to. It also has what high level codes correspond to what low level codes we use here (and our lumping/splitting schema).

studyextent=readOGR(dsn = "cloudintersect_AKALB.shp", layer = "cloudintersect_AKALB") #Loads a layer with the Cropping mask for a) study area and b) clouds that were in the imgagery, so we don't sample 'truth' points from there.

#Random sample pointn points for each LC (sampleStratified) from the brick so we have pointn instances of each LC, and their corresponding enviromental covariates. 
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

#Here we join our LC descriptions to numbers. Important for tracking things through the code.
model_data=randompts_data@data
model_data$LC=join(randompts_data@data,LCCODEX,by="LC",type="left")$LCrecode
model_data=model_data[order(model_data$LC),] ## This does not always work due to memory limitations. 

#Train a model using random forests.
Sys.time()
mark=Sys.time()
modelRF = train(as.factor(LC) ~ B2 + B3 + B4 + B5 + B6 + B7 + B10 + B11 + SOLAR + DEM + SLP + IFSAR + S2 + S3 + S4 + S5 + S6 + S7 + S10 + S11, method = "rf", data=subset(model_data[complete.cases(model_data),], (LC != 20)), tuneLength = 6) #Note lack of B1 (Aerosol) B8 (Panchro) and B9 (Cirrus). data statement is to remove NAs and invalid LCs (recoded as 20s), as well as any pts that happened off the raster.
mark2=Sys.time()
modelfittime=mark2-mark
modelRF

#Predict a new raster based on the final RF model. I use Cluster to speed things along considerably.
Sys.time()
mark=Sys.time()
beginCluster()
preds_rf <- clusterR(envcov, raster::predict, args = list(model = modelRF))
endCluster()
mark2=Sys.time()
predicttime=mark2-mark

#Write our results to disk.
writeRaster(preds_rf,"LC_OUTPUT_30M_FS_9-22.tif")
