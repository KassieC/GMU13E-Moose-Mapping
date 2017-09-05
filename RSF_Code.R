library(plyr)
library(sp)
library(raster)
library(rgdal)
wd="" ##Directory with the data in it.

setwd(wd)
LOC_DF=read.table("WatMooPts-1267.txt",header=TRUE,sep="\t")
masterid=read.table("masteridlist.txt",header=TRUE,sep="\t")

LOC_DF$jday=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$yday+1
LOC_DF$hour=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$hour+(as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$min/60)+(as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$sec/(60*60))
LOC_DF$year=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$year+1900
LOC_DF$gmt=as.POSIXct(strptime(paste(LOC_DF$Acquisition_Time),"%Y.%m.%d %H:%M:%S"), tz="GMT")
LOC_DF$AnNum=gsub(",","",LOC_DF$AnNum)
LOC_DF$sex=join(LOC_DF,masterid,by="AnNum",type="left",match="all")$Sex
LOC_DF$Cohort=as.factor(LOC_DF$Cohort)

## Calendar
## if jday >= 336 OR <= 90 season is winter ##Bracket Dec1 to Mar31
## if jday >= 91 AND <= 129 is spring ##Bracket Apri1 to May9
## if jday >= 130 AND <= 166 is calving ##Bracket May 10 to June 15
## if jday >= 167 AND <= 243 is summer ##Bracket June 16 Aug 31
## if jday >= 244 AND <= 304 is fall ##Bracket Sept 1 to Oct 31
## if jday >= 305 AND <= 335 is post-rut ##Bracket Nov 1 to Nov 30

season=c(rep("latewinter",91),rep("spring",39),rep("calving",37),rep("summer",77),rep("fall",61),rep("postrut",30),rep("winter",31)) #This could probably be handled better with an ifelse. Calendar is 366 days long (leap year).
cal=as.data.frame(x = cbind(1:366,season)) 
cal$season=as.factor(cal$season)
names(cal)=c("jday","season")
LOC_DF$eseason=join(LOC_DF,cal,by="jday",type="left",match="all")$season #Extended season. For housekeeping purposes.
LOC_DF$pseason=as.factor(ifelse(LOC_DF$eseason == "latewinter",paste("winter",sep=""),(paste(LOC_DF$eseason,sep="")))) #Pooled season. Pools season across all years.
LOC_DF$yseason=as.factor(ifelse(LOC_DF$eseason == "latewinter",(paste(LOC_DF$year-1,"winter",sep="")),(paste(LOC_DF$year,LOC_DF$eseason,sep="")))) #Year-season, wrapping winters across the year boundary, and then storing it by-year

# LOC_DF$yearseasonid=as.factor(paste(LOC_DF$AnNum,LOC_DF$yseason,sep="_")) # This is only needed if we're doing KDEs by year-season-annum

#Mark these as used points.
LOC_DF$used=rep(1,nrow(LOC_DF))

#Write drop table
DTable=as.data.frame(cbind(ifelse(test = (summary(as.factor(LOC_DF$AnNum))<=200),1,0),names(summary(as.factor(LOC_DF$AnNum)))))
row.names(DTable)=c()
names(DTable)=c("drop","AnNum")
LOC_DF$drop=join(LOC_DF,DTable,by="AnNum",type="left",match="all")$drop
#Dropping any animals where n(total) < 200
LOC_DF=subset(LOC_DF,(drop==0))
LOC_DF=droplevels(LOC_DF)


#Make this into a spatial df
LOC_SPDF<-SpatialPointsDataFrame(coords = cbind(LOC_DF$GPS_Longitude,LOC_DF$GPS_Latitude), data=LOC_DF, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
LOC_SPDF$AnNum=as.factor(LOC_SPDF$AnNum)
LOC_SPDF$yseason=as.factor(LOC_SPDF$yseason)
LOC_SPDF$pseason=as.factor(LOC_SPDF$pseason)


#####Loading Rasters###############

rasterOptions(tmpdir=paste(wd,"/tempdir",sep = "")) #this is important to dealing with raster grids... because otherwise they are going to vanish in the default tempfile sometime during your analysis. That's a bad thing(tm)
dem=raster("DEM_RS_Scaled.tif")
dem2=calc(dem, fun=function(x){x^2})
TRI=raster("TRI_Scaled.tif")
LC_O=raster("LC_OUTPUT_30M_FS_9-22b.tif")

SLP=raster("SLP_RS.tif")
# aspect=raster("ASP_RS.tif")
# transformed_aspect=calc(x = aspect, fun = sin)
# writeRaster(transformed_aspect,"ASP_SIN_RS.tif")
# slope_aspect_interaction = slope * transformed_aspect
# writeRaster(slope_aspect_interaction,"ASP_INT_RS.tif")
ASP=raster("ASP_SIN_RS.tif")
ASPINT=raster("ASP_INT_RS.tif")

LCRECODE=read.table("recodelist.txt",header=TRUE,sep="\t")
LC = reclassify(LC_O,rcl = cbind(LCRECODE$LC_II_recode,LCRECODE$LC_RSF)[!duplicated(cbind(LCRECODE$LC_II_recode,LCRECODE$LC_RSF)),])

LC = as.factor(LC)
temp = levels(LC)[[1]]
temp$code=c("OPN", "NDL", "H2O", "BDL", "LSC", "TSC", "MIX", "BAR")
levels(LC)=temp
rm(temp)
LC_brick=stack(LC)
for (i in 1:length(levels(LC)[[1]][,1]))
  {
  LC_temp=calc(LC,fun=function(x){ifelse(x==levels(LC)[[1]][i,1],1,0)})
  LC_brick[[i]]=LC_temp
  names(LC_brick)[i]=levels(LC)[[1]][i,2]
  }

bricked=stack(dem,dem2,TRI,LC_brick,LC_O,SLP,ASP,ASPINT)
names(bricked)=c("dem","dem2","TRI",names(LC_brick),"LC_O","SLP","ASP","ASPINT")


i=1
LOC_SPDF=spTransform(LOC_SPDF,bricked@crs) #Re-project to AK Albers
#Why AK Albers? Because the final product is a map where we're going to look at the habitat value of different parcels. In this context, an equal area project seems best. Defining now is important since all of this is getting passed to JAGS for predicting, so we need our matrix set up to meet our end needs before its passed.

availarea=readOGR(dsn = "1267AllBuff/1267AllBuff.shp", layer = "1267AllBuff")
availarea=spTransform(availarea,bricked@crs)

#####Sampling Available points############
# This is an odd bit of code here. It comes from the fact that if-else statements evaluate both the then and the else, regardless of what is going to be logically chosen, so if you have an invalid statement in the 'else' statement, the whole thing falls apart because package::sp has horrible error handling and would rather see the world burn than cleanly fail. See http://stackoverflow.com/questions/16275149/does-ifelse-really-calculate-both-of-its-vectors-every-time-is-it-slow Normally I'd just do nested whiles for j seasons and i individuals and sample random points from homerange i when length(seasonXindividual) > 0, but R isn't going for that because when length zero, the spsample command throws an exception that brings things to a screaming halt. Instead, I've created a bit of code that won't throw exceptions. In a nutshell, it creates an index object of animals and seasons that do exist. In a second while loop (which should be a for, but that doesn't work for ??? reasons) I move along the rows of the index object to retrieve the i and j for animals that will work without throwing an exception.
index=data.frame(i=integer(),j=integer())
k=1
for (j in 1:length(levels(LOC_SPDF@data$pseason)))
{
  for(i in 1:length(levels(LOC_SPDF@data$AnNum)))
  {
    dummy=eval(nrow(subset(LOC_SPDF,(LOC_SPDF@data$AnNum==levels(LOC_SPDF@data$AnNum)[i]) & (LOC_SPDF@data$pseason==levels(LOC_SPDF@data$pseason)[j])))==0)
    if (dummy==TRUE)
    {cat(paste(levels(LOC_SPDF@data$AnNum)[i],levels(LOC_SPDF@data$pseason)[j],"does not exist.\n",sep=" "))}
    else
    {index[k,]=c(i,j);k=k+1} 
  }
}


# A few brief notes. This is currently set up to select random points from the pooled 12.67km buffered area around all points (2nd order selection - though this is fuzzy). It does this for inf_factor times the number of used points in a given season (so if pooled spring has 300 points, it samples 900 avail points from the availarea poly). That number can be easily changed to any number of avail. points where the dumpspdf object is being written in the first line. dumpspdf, as the name implies, is a temporary dump object that gets bound to the temporary AVAIL_SPDF. 12.67 is the radius of a circle that leads to a 200 mile^2 area (approximately, anyhow), which is about the HR found by Ballard in GMU 13B. This is a generalization, but since we're looking at home-range selection, it should suffice. On top of this, I've eliminated all points within 12.67 km of the boundary of the rasters, because these would have potentail avail points outside the WSA rasters and would have NA data. This is why sample size for some of the individuals is so low, because they've had most of their points dropped. Luckily, we start with an absurd # of individuals so the loss of a few isn't a big deal. 


inf_factor = 8 #Number of avail-per-used

AVAIL_SPDF=LOC_SPDF[1,c("AnNum", "sex", "pseason", "used")]
AVAIL_SPDF=AVAIL_SPDF[-1,]


for (i in 1:nrow(LOC_SPDF))
  {dumpspdf = spsample(buffer(LOC_SPDF[i,],width=12670), n = inf_factor, type = "random", iter=10)
  dumpdata=as.data.frame(cbind(rep((levels(LOC_SPDF@data$AnNum)[LOC_SPDF[i,]$AnNum]),inf_factor),rep((levels(LOC_SPDF@data$sex)[LOC_SPDF[i,]$sex]),inf_factor),rep(levels(LOC_SPDF@data$pseason)[LOC_SPDF[i,]$pseason],inf_factor),rep(0,inf_factor)))
  names(dumpdata)=c("AnNum","sex","pseason","used")  
  dumpspdf=SpatialPointsDataFrame(dumpspdf,data=dumpdata)
  AVAIL_SPDF=rbind(AVAIL_SPDF,dumpspdf)
  }
rm(dumpspdf,dumpdata,i)

# In order to add the sex information, I need to unload adehabitat* (and re-load plyr) to unmask plyr's ID function. Frustrating. 
# detach(package:adehabitatLT)
# detach(package:adehabitatMA)
# detach(package:adehabitatHR)
# detach(package:plyr)
# library(plyr)

temp=as.data.frame(x= c(as.character(AVAIL_SPDF$AnNum)))
names(temp)=c("AnNum")
AVAIL_SPDF$sex=join(temp,masterid,by="AnNum",type="left",match="all")$Sex
rm(temp) 

#This is just some error checking code to make sure I'm actually generating my avail correctly on a per-animal, per-pooledseason basis. It should be commented out when not in use.
# temp=cbind(c(),c())
# for(i in 1:length(levels(LOC_SPDF$AnNum)))
# {temp=rbind(temp,c(length(subset(LOC_SPDF, (AnNum == levels(LOC_SPDF$AnNum)[i])))/length(subset(AVAIL_SPDF, (AnNum == levels(LOC_SPDF$AnNum)[i]))), length(subset(LOC_SPDF, (AnNum == levels(LOC_SPDF$AnNum)[i]))), length(subset(AVAIL_SPDF, (AnNum == levels(LOC_SPDF$AnNum)[i]))), levels(LOC_SPDF$AnNum)[i]))}

#####Merge db and extract.##############
ANALYSIS_SPDF=rbind(LOC_SPDF[,c("AnNum","sex","pseason","used")],AVAIL_SPDF) #Merges only the AnNum, sex, pseason, and used colums ATM. In the future, we may wish to include whether an animal was parturant or not for calving model (YIKES!)
ANALYSIS_SPDF@data$used=as.numeric(ANALYSIS_SPDF@data$used)

ANALYSIS_SPDF@data=cbind(ANALYSIS_SPDF@data, extract(bricked,ANALYSIS_SPDF))

cor.test(ANALYSIS_SPDF@data$TRI,ANALYSIS_SPDF@data$SLP,method="pearson") #Slope and TRI correlated.
cor.test(ANALYSIS_SPDF@data$ASP,ANALYSIS_SPDF@data$SLP,method="pearson") #Slope and ASP not correlated.
cor.test(ANALYSIS_SPDF@data$ASP,ANALYSIS_SPDF@data$TRI,method="pearson") #TRI and ASP not correlated.
cor.test(ANALYSIS_SPDF@data$ASPINT,ANALYSIS_SPDF@data$TRI,method="pearson") #TRI and ASP not correlated.

m1 = glm(used ~ TRI, data= ANALYSIS_SPDF@data, family = "binomial")
m2 = glm(used ~ SLP, data= ANALYSIS_SPDF@data, family = "binomial")
m3 = glm(used ~ TRI + ASP, data= ANALYSIS_SPDF@data, family = "binomial")
m4 = glm(used ~ TRI + ASP + TRI*ASP, data= ANALYSIS_SPDF@data, family = "binomial")
m5 = glm(used ~ SLP + ASP, data= ANALYSIS_SPDF@data, family = "binomial")
m6 = glm(used ~ SLP + ASP + ASPINT, data= ANALYSIS_SPDF@data, family = "binomial")
AIC(m1,m2,m3,m4,m5,m6) #Model 4 (TRI, ASP, and an interaction term) is best supported. TRI outperforms slope in all models. ASP is in top models for both slope and TRI.

ASPINT2 = ASP * TRI

m1=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1,data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)), family = "binomial")
m2=glmer(used~dem+dem2+SLP+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1,data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)), family = "binomial")
m3=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+ASP+TRI*ASP+(1|AnNum)-1,data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)), family = "binomial")
m4=glmer(used~dem+dem2+SLP+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+ASP+SLP*ASP+(1|AnNum)-1,data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)), family = "binomial")
m5=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+ASP+(1|AnNum)-1,data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)), family = "binomial")
AIC(m1,m2,m3,m4,m5)

#Models.
library(lme4)
library(lmerTest)
library(coefplot2)
library(ResourceSelection)
library(snow)
library(boot)
Sys.time()



kfolds=5
ModSeason="winter"
ModSex="M"
gen_form="used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1" #This is where we're storing the glmer formula. Will paste it elsewhere to save on linespace and to allow easy editing of all models at the same time. call as "as.formula(gen_form)"
#gen_form="used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+ASP+TRI*ASP+(1|AnNum)-1"
#model=glmer(as.formula(gen_form),data=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason=="calving"|pseason=="summer")), family="binomial",nAGQ=6)

temp_dat_obj=subset(ANALYSIS_SPDF@data,(sex==ModSex)&(pseason==ModSeason)) #creating subset data object for the model in question.
withhold_key=as.data.frame(cbind(levels(temp_dat_obj$AnNum), sample(1:kfolds, length(levels(temp_dat_obj$AnNum)), replace=TRUE))) ## Note k-folds is conditioned on animal (annum), using 5 folds.
names(withhold_key)=c("AnNum","Fold_N")
temp_dat_obj$Fold_N=join(temp_dat_obj,withhold_key,by="AnNum",type="left",match="all")$Fold_N ## just joining up the key to datasubset so everyone has their fold number
rank_cor_obj=c()
for(i in 1:kfolds)
  {
  trained=subset(temp_dat_obj,(Fold_N!=i))  #Training data = everyone but fold i
  witheld=subset(temp_dat_obj,(Fold_N==i))  #Withheld data = only fold i
  predicted_model=glmer(as.formula(gen_form),data=trained, family="binomial")
  beginCluster()
  predicted_rast <- clusterR(bricked, raster::predict, args = list(model = predicted_model, re.form = NA)) #Predict raster based on model with no randomeffects structure.
  endCluster()
  predicted_witheld=predict(predicted_model, newdata=subset(witheld,(used==1)),re.form=NA) #Predict RSF scores to the withheld used points. This is used later...
  predicted_avail_fit=extract(predicted_rast,availarea) #extract from the totality of the study area.
  ###Following code care of Katie Christie/Mark Boyce, edited to play nice with us.####
  ###Biggest changes are generalization to work with looping arch. as well as storing outputs to df for automation purposes.
  names(predicted_avail_fit)="fit"
  pred.avail.fit<-na.omit(predicted_avail_fit$fit)
  #Make rsf bins with equal numbers of random locations in each bin
  avail.rsf.bins<-cut(pred.avail.fit, quantile(pred.avail.fit, seq(0, 1, len = 11)), include.lowest = TRUE) 
  avail.rsf.bins<-table(avail.rsf.bins)
  avail.rsf.bins
  v<-dimnames(avail.rsf.bins)
  v <- data.frame(v, stringsAsFactors = FALSE)
  y1<-v[1,1]
  y2<-v[3,1]
  y3<-v[5,1]
  y4<-v[7,1]
  y5<-v[9,1]
  y6<-v[10,1]
  ynew<-c(y1,y2,y3,y4,y5,y6)
  ynew<-gsub("\\[", "", ynew)
  ynew<-gsub("\\]", "", ynew)
  ynew<-gsub("\\(", "", ynew)
  ynew<-gsub("  ", "", ynew)
  ynew<-as.numeric(unlist(strsplit(ynew,split="[,]")))
  ynew<-ynew[-10]
  avail.rsf.bins<-ynew
  
  rsf.with<-predicted_witheld[1:length(predicted_witheld)]
  hist(rsf.with)
  rsf.bins<-cut(rsf.with,breaks=avail.rsf.bins)
  rsf.bins<-table(rsf.bins)
  rsf.bins<-as.numeric(rsf.bins)
  binrank<-seq(1,10,length.out=10)
  cor_out=cor.test(rsf.bins,binrank,method="spearman")
  rank_cor_obj=rbind(rank_cor_obj,c(rsf.bins,cor_out$p.value,cor_out$estimate))
  cat(paste("Fold",i,"of",kfolds))
  cat(Sys.time())
  }

rank_cor_overall=cor.test(c(rank_cor_obj[1,1:10],rank_cor_obj[2,1:10],rank_cor_obj[3,1:10],rank_cor_obj[4,1:10],rank_cor_obj[5,1:10]),rep(1:10,5),method="spearman")

Winter_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="winter")), family="binomial", nAGQ=3)
Sys.time()
Spring_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="spring" | pseason=="calving")), family="binomial", nAGQ=3)
Sys.time()
# Calving_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="calving")), family="binomial", nAGQ=3)
# Sys.time()
Summer_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="summer")), family="binomial", nAGQ=3)
Sys.time()
Fall_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="postrut" | pseason=="fall")), family="binomial", nAGQ=3)
Sys.time()
# Postrut_model=glmer(used~dem+dem2+TRI+OPN+NDL+H2O+BDL+LSC+TSC+MIX+BAR+(1|AnNum)-1 ,data=subset(ANALYSIS_SPDF@data,(sex=="M")&(pseason=="postrut")), family="binomial", nAGQ=3)

par(mfrow=c(2,3))

####TRI Residual plot
plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="winter"))$TRI, residuals(Winter_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="winter"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="spring"))$TRI, residuals(Spring_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="spring"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="calving"))$TRI, residuals(Calving_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="calving"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="summer"))$TRI, residuals(Summer_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="summer"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="fall"))$TRI, residuals(Fall_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="fall"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="postrut"))$TRI, residuals(Postrut_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="postrut"))$used+1])

####Elevation Residual plot
plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="winter"))$dem, residuals(Winter_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="winter"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="spring"))$dem, residuals(Spring_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="spring"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="calving"))$dem, residuals(Calving_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="calving"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="summer"))$dem, residuals(Summer_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="summer"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="fall"))$dem, residuals(Fall_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="fall"))$used+1])

plot(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="postrut"))$dem, residuals(Postrut_model), col=c("Blue","Red")[subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="postrut"))$used+1])


beginCluster()
Fall_Pred <- clusterR(bricked, raster::predict, args = list(model = Fall_model, re.form = NA))
Winter_Pred <- clusterR(bricked, raster::predict, args = list(model = Winter_model, re.form = NA))
Postrut_Pred <- clusterR(bricked, raster::predict, args = list(model = Postrut_model, re.form = NA))
endCluster()
Winter_Pred = calc(Winter_Pred, fun=inv.logit)
Postrut_Pred = calc(Postrut_Pred, fun=inv.logit)
Fall_Pred= calc(Fall_Pred, fun=inv.logit)

par(mfrow=c(1,2))
plot(Winter_Pred)
plot(Postrut_Pred)
par(mfrow=c(1,1))

hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="winter"))$used, fitted(Winter_model), 10)
hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="spring"))$used, fitted(Spring_model), 10)
hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="calving"))$used, fitted(Calving_model), 10)
hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="summer"))$used, fitted(Summer_model), 10)
hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="fall"))$used, fitted(Fall_model), 10)
hoslem.test(subset(ANALYSIS_SPDF@data,(sex=="F")&(pseason=="postrut"))$used, fitted(Postrut_model), 10)
