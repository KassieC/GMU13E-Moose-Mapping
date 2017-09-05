library(plyr)
library(sp)
library(raster)
library(rgdal)
wd="" #Directory with data

setwd(wd)
LOC_DF=read.csv("WatanaMoosePts.txt",header=TRUE)
masterid=read.table("masteridlist.txt",header=TRUE,sep="\t")
AKALBCRS=raster("finalBrick.tif")@crs


LOC_DF$jday=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$yday+1
LOC_DF$hour=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$hour+(as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$min/60)+(as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$sec/(60*60))
LOC_DF$year=as.POSIXlt(as.POSIXct(LOC_DF$Acquisition_Time, tz="GMT", format="%Y.%m.%d %H:%M:%S"), tz="US/Alaska")$year+1900
LOC_DF$gmt=as.POSIXct(strptime(paste(LOC_DF$Acquisition_Time),"%Y.%m.%d %H:%M:%S"), tz="GMT")
LOC_DF$AnNum=gsub(",","",LOC_DF$AnNum)
LOC_DF$sex=join(LOC_DF,masterid,by="AnNum",type="left",match="all")$Sex
#LOC_DF$Cohort=as.factor(LOC_DF$Cohort)

## Calendar
## if jday >= 336 OR <= 90 season is winter ##Bracket Dec1 to Mar31
## if jday >= 91 AND <= 129 is spring ##Bracket Apri1 to May9
## if jday >= 130 AND <= 166 is calving ##Bracket May 10 to June 15
## if jday >= 167 AND <= 243 is summer ##Bracket June 16 Aug 31
## if jday >= 244 AND <= 304 is fall ##Bracket Sept 1 to Oct 31
## if jday >= 305 AND <= 335 is post-rut ##Bracket Nov 1 to Nov 30

season=c(rep("latewinter",91),rep("spring",39),rep("calving",37),rep("summer",77),rep("fall",61),rep("postrut",30),rep("winter",31)) #This could probably be handled better with an ifelse. Calendar is 366 days long (leap year).
cal=as.data.frame(x = cbind(1:366,season)) 
names(cal)=c("jday","season")
cal$season=as.factor(cal$season)
cal$jday=as.numeric(1:366)
LOC_DF$eseason=join(LOC_DF,cal,by="jday",type="left",match="first")$season #Extended season. For housekeeping purposes.
LOC_DF$pseason=as.factor(ifelse(LOC_DF$eseason == "latewinter",paste("winter",sep=""),(paste(LOC_DF$eseason,sep="")))) #Pooled season. Pools season across all years.
LOC_DF$yseason=as.factor(ifelse(LOC_DF$eseason == "latewinter",(paste(LOC_DF$year-1,"winter",sep="")),(paste(LOC_DF$year,LOC_DF$eseason,sep="")))) #Year-season, wrapping winters across the year boundary, and then storing it by-year

LOC_DF$yearseasonid=as.factor(paste(LOC_DF$AnNum,LOC_DF$yseason,sep="_")) 

LOC_DF$yearid=as.factor(paste(LOC_DF$AnNum,LOC_DF$year,sep="_")) 

part_table=read.delim("parturition.txt",header = TRUE, sep="\t")
part_table$yearseasonid=as.factor(paste(part_table$AnNum,"_",part_table$Year,"calving",sep=""))
part_table$yearid=as.factor(paste(part_table$AnNum,"_",part_table$Year,sep=""))

LOC_DF$partur=join(LOC_DF,part_table,by="yearseasonid",type="left",match="first")$Parturition
LOC_DF$partur2=join(LOC_DF,part_table,by="yearid",type="left",match="first")$Parturition


#What are we doing here? Origionally our calendar contained 6 seasons (and winter split), but we need to relump everything to fix some issues we've had. Unfortunately, we can't alter the calendar directly (legacy reasons), and in any event, the calendar can't help with parturition. So, for purposes here, we're going to chop up the dataset and re-classify everything into new seasons manually.
temp_df1=subset(LOC_DF,(pseason=="winter")|(pseason=="summer")) #These are the points we don't want to re-classify

temp_df2=subset(LOC_DF,(pseason=="fall")|(pseason=="postrut")) #These are the points we're lumping into "autumn"
temp_df2$pseason=rep("autumn",nrow(temp_df2))

temp_df3=subset(LOC_DF,((pseason=="spring")|(pseason=="calving"))&(sex=="M")) #These are the male points we're calling "greenup"
temp_df3$pseason=rep("greenup",nrow(temp_df3))

temp_df4=subset(LOC_DF,(pseason=="spring")&(sex=="F")) #All female spring points are re-classed "greenup"
temp_df4$pseason=rep("greenup",nrow(temp_df4))            

temp_df5=subset(LOC_DF,(pseason=="calving")&(sex=="F")&((partur==0)|(is.na(partur)==TRUE))) #These are the calf-less female points we're calling "greenup"
temp_df5$pseason=rep("greenup",nrow(temp_df5))  

temp_df6=subset(LOC_DF,(pseason=="calving")&(sex=="F")&(partur==1)) #These are the parturant female points we're calling "calving" They start out as calving, so don't need over-written

LOC_DF=rbind(temp_df1,temp_df2,temp_df3,temp_df4,temp_df5,temp_df6) #now we rebind the whole lot into an object
rm(temp_df1,temp_df2,temp_df3,temp_df4,temp_df5,temp_df6) #and erase the temporary files.
LOC_DF=droplevels(LOC_DF) #Trim un-used levels.


DTable=subset(aggregate(jday~AnNum, data=LOC_SPDF@data, NROW),(jday>=365*3))
row.names(DTable)=c()
names(DTable)=c("AnNum","DayCount")
DTable$drop=c(rep(0,nrow(DTable)))
LOC_DF$drop=join(LOC_DF,DTable,by="AnNum",type="left",match="all")$drop
#Dropping any animals where there is less than 365 days worth of data
#LOC_DF=subset(LOC_DF,(drop==0))
LOC_DF=droplevels(LOC_DF)
#Remove 2016 calving and summer as they're outside the study period.
LOC_DF=subset(LOC_DF,(yseason!="2016calving")&(yseason!="2016summer")) #Note this doesn't bring anyone new below 365*3

LOC_SPDF<-SpatialPointsDataFrame(coords = cbind(LOC_DF$GPS_Longitude,LOC_DF$GPS_Latitude), data=LOC_DF, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
LOC_SPDF$AnNum=as.factor(LOC_SPDF$AnNum)
LOC_SPDF$yseason=as.factor(LOC_SPDF$yseason)
LOC_SPDF$pseason=as.factor(LOC_SPDF$pseason)
LOC_SPDF$year_fac=as.factor(LOC_SPDF$year)

LOC_SPDF=spTransform(LOC_SPDF,AKALBCRS)
LOC_SPDF@data$Frac_Year=LOC_SPDF@data$year+(LOC_SPDF@data$jday/366)+((LOC_SPDF@data$hour/24)/366)
LOC_SPDF@data$Frac_Day=LOC_SPDF@data$jday+(LOC_SPDF@data$hour/24)

LOC_SPDF=LOC_SPDF[order(LOC_SPDF@data$AnNum,LOC_SPDF@data$year,LOC_SPDF@data$Frac_Day),] #Sort the df so it is ordered by animals and then by time 
LOC_SPDF$velocity=as.numeric(rep(NA,nrow(LOC_SPDF)))


vel_SPDF=LOC_SPDF[1,]
vel_SPDF=vel_SPDF[-1,]


for (i in 1:length(levels(LOC_SPDF$AnNum)))
  {
  temp_SPDF=subset(LOC_SPDF,(AnNum==levels(LOC_SPDF$AnNum)[i]))
  temp_SPDF@data=droplevels(temp_SPDF@data)
  for (j in 1:length(levels(temp_SPDF$year_fac)))
    {
    temp_SPDF_2=subset(temp_SPDF,(year_fac==levels(temp_SPDF$year_fac)[j]))
    for (k in 1:(nrow(temp_SPDF_2)-1))
      {
      temp_SPDF_2$velocity[k+1]=(spDistsN1(temp_SPDF_2[k,],temp_SPDF_2[k+1,])/((temp_SPDF_2$Frac_Day[k+1]-temp_SPDF_2$Frac_Day[k])*24))
      }
    vel_SPDF=rbind(vel_SPDF,temp_SPDF_2)
    }
  }
plot(aggregate(velocity~jday, data=subset(vel_SPDF@data,(sex=="M")), FUN=mean),type="l",col="red",log="y")
points(aggregate(velocity~jday, data=subset(vel_SPDF@data,(sex=="F")), FUN=mean),type="l",col="black")

#Create seasonal plots.
weightedavgvel=as.data.frame(cbind(c(1:366),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366))))
names(weightedavgvel)=c("jday","cow.vel","bull.vel","cow.Q1","bull.Q1","cow.Q3","bull.Q3")

#Add days upto -10 and 376 to the calendar by wrapping around the newyears. Allows us to do the moving window.
temp=vel_SPDF
temp2=subset(vel_SPDF,(jday >= 366-10))
temp2$jday=temp2$jday - 366
temp=rbind(temp,temp2)
rm(temp2)
temp2=subset(vel_SPDF,(jday <= 10))
temp2$jday=temp2$jday + 366
temp=rbind(temp,temp2)
rm(temp2)


temp=subset(temp,(velocity<=10000))
for (i in 1:length(levels(vel_SPDF$sex)))
  {for (j in 1:nrow(weightedavgvel))
    {
      weightedavgvel[j,i+1]=median(subset(temp,(sex == levels(temp$sex)[i])&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, na.rm = TRUE)
      weightedavgvel[j,i+3]=quantile(subset(temp,(sex == levels(temp$sex)[i])&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.25, na.rm = TRUE)
      weightedavgvel[j,i+5]=quantile(subset(temp,(sex == levels(temp$sex)[i])&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.75, na.rm = TRUE)
    }  
  }



greenupstart=min(subset(temp,(sex=="F")&(pseason=="greenup"))$jday)
greenupend=max(subset(temp,(sex=="F")&(pseason=="greenup"))$jday)
summerend=max(subset(temp,(sex=="F")&(pseason=="summer"))$jday)
calvingstart=min(subset(temp,(sex=="F")&(pseason=="calving"))$jday)

weighted_f_calv=as.data.frame(cbind(c(1:366),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366)),c(rep(NA,366))))
names(weighted_f_calv)=c("jday","greenup.vel","calv.vel","greenup.Q1","calv.Q1","greenup.Q3","calv.Q3")


####This code just does Greenup+Calving and returns summer to normal. Commented out because parturient cows moved less in summer.
# temp=subset(temp,(sex=="F")&((jday>=greenupstart)&(jday<=summerend)))
# temp@data=droplevels(temp@data)
# for (i in 1:length(levels(temp$pseason)))
# {
# for (j in greenupstart:summerend)
# {
#   weighted_f_calv[j,i+1]=median(subset(temp,(pseason == levels(temp$pseason)[i]&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, na.rm = TRUE)
#   weighted_f_calv[j,i+3]=quantile(subset(temp,(pseason == levels(temp$pseason)[i])&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.25, na.rm = TRUE) 
#   weighted_f_calv[j,i+5]=quantile(subset(temp,(pseason == levels(temp$pseason)[i])&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.75, na.rm = TRUE)
#   }
# }

temp=subset(temp,(velocity<=10000))
temp=subset(temp,(sex=="F")&((jday>=greenupstart)&(jday<=summerend)))
temp@data=droplevels(temp@data)
i=1
for (j in greenupstart:summerend)
{
  weighted_f_calv[j,i+1]=median(subset(temp,(partur2!=1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, na.rm = TRUE)
  weighted_f_calv[j,i+3]=quantile(subset(temp,(partur2!=1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.25, na.rm = TRUE) 
  weighted_f_calv[j,i+5]=quantile(subset(temp,(partur2!=1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.75, na.rm = TRUE)
}
i=2
for (j in greenupstart:summerend)
{
  weighted_f_calv[j,i+1]=median(subset(temp,(partur2==1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, na.rm = TRUE)
  weighted_f_calv[j,i+3]=quantile(subset(temp,(partur2==1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.25, na.rm = TRUE) 
  weighted_f_calv[j,i+5]=quantile(subset(temp,(partur2==1)&(((jday <= j+3)&(jday >= j))|((jday >= (j-3))&(jday <= j))))$velocity, probs=.75, na.rm = TRUE)
}





# weighted_f_calv$calv.vel[weighted_f_calv$jday < calvingstart ] = NA #Doesn't do anything but such a clean bit of code I can't bring myself to delete it.
# weighted_f_calv$calv.05[weighted_f_calv$jday < calvingstart ] = NA 
# weighted_f_calv$calv.95[weighted_f_calv$jday < calvingstart ] = NA

weighted_f_calv=subset(weighted_f_calv,(jday>=calvingstart))

maxcowplot=max(c(weightedavgvel$cow.Q3,weighted_f_calv$greenup.Q3),na.rm=TRUE)

par(mfrow=c(2,1))
plot(weightedavgvel$jday,weightedavgvel$cow.vel,type="l",ylim=c(0,maxcowplot),xlab="",ylab="Net movement since previous relocation (m/hr)",main="Cow moose seasonal movement rate")
rect(0,0,91,maxcowplot,col="lightblue", border=NA)
rect(91,0,130,maxcowplot,col="lightgreen", border=NA)
rect(130,0,167,maxcowplot,col="lightcoral", border=NA)
rect(167,0,244,maxcowplot,col="green4", border=NA)
rect(244,0,305,maxcowplot,col="yellow", border=NA)
rect(305,0,336,maxcowplot,col="yellow", border=NA)
rect(336,0,366,maxcowplot,col="lightblue", border=NA)

points(weightedavgvel$jday,weightedavgvel$cow.vel,type="l")
points(weightedavgvel$jday,weightedavgvel$cow.Q1,type="l",lty=2)
points(weightedavgvel$jday,weightedavgvel$cow.Q3,type="l",lty=2)
#rect(91,0,130,max(weightedavgvel$cow.Q3),col="lightgreen", border=NA)
rect(130,0,167,maxcowplot,col="lightcoral", border=NA)
rect(167,0,244,maxcowplot,col="green4", border=NA)


points(weighted_f_calv$jday,weighted_f_calv$greenup.vel,type="l")
points(weighted_f_calv$jday,weighted_f_calv$greenup.Q1,type="l",lty=2)
points(weighted_f_calv$jday,weighted_f_calv$greenup.Q3,type="l",lty=2)

points(weighted_f_calv$jday,weighted_f_calv$calv.vel,type="l",lwd=2,col="firebrick")
points(weighted_f_calv$jday,weighted_f_calv$calv.Q1,type="l",lty=2, lwd=2,col="firebrick")
points(weighted_f_calv$jday,weighted_f_calv$calv.Q3,type="l",lty=2, lwd=2,col="firebrick")

mtext(c("Winter","Green-up","Calving","Summer","Autumn","Winter"),cex=1.25,at=c(55,110,150,200,300,350),side=1,line=2)
mtext(c("or Green-Up"),cex=1.25,at=c(150),side=1,line=3)

plot(weightedavgvel$jday,weightedavgvel$bull.vel,type="l",ylim=c(0,max(weightedavgvel$bull.Q3)),xlab="Julian Date",ylab="Net movement since previous relocation (m/hr)",main="Bull moose seasonal movement rate")
rect(0,0,91,max(weightedavgvel$bull.Q3),col="lightblue", border=NA)
rect(91,0,130,max(weightedavgvel$bull.Q3),col="lightgreen", border=NA)
rect(130,0,167,max(weightedavgvel$bull.Q3),col="lightgreen", border=NA)
rect(167,0,244,max(weightedavgvel$bull.Q3),col="green4", border=NA)
rect(244,0,305,max(weightedavgvel$bull.Q3),col="yellow", border=NA)
rect(305,0,336,max(weightedavgvel$bull.Q3),col="yellow", border=NA)
rect(336,0,366,max(weightedavgvel$bull.Q3),col="lightblue", border=NA)

points(weightedavgvel$jday,weightedavgvel$bull.vel,type="l")
points(weightedavgvel$jday,weightedavgvel$bull.Q1,type="l",lty=2) 
points(weightedavgvel$jday,weightedavgvel$bull.Q3,type="l",lty=2) 

mtext(c("Winter","Green-up","Summer","Autumn","Winter"),cex=1.25,at=c(55,120,200,300,350),side=1,line=2)


par(mfrow=c(1,1))

vel_SPDF@data$pseason=factor(vel_SPDF@data$pseason, levels=c("winter", "greenup", "calving", "summer", "autumn")) #reorder the seasons as factors.
boxplot(velocity~sex*pseason,data=subset(vel_SPDF,(velocity<=10000)), col=(c("lightblue","cyan4","lightgreen","springgreen4","lightcoral","lightcoral","olivedrab3","olivedrab4","yellow","gold")),outline=FALSE,names=c("Cow","Bull","Cow","Bull","Cow","Bull","Cow","Bull","Cow","Bull"),ylab="Net movement since previous relocation (m/hr)")
text(6,50,labels="NA",cex=2)
mtext(c("Winter","Green-up","Calving","Summer","Autumn"),cex=1.5,at=(c(1:5)*2)-0.5,side=1,line=3) 



#Dropping any animals where there is less than 365 days worth of data
LOC_DF=subset(LOC_DF,(drop==0))


library(ade4)
library(adehabitatHR)
library(adehabitatHS)
library(adehabitatLT)
library(adehabitatMA)

kud_SPDF=SpatialPointsDataFrame(coords = cbind(LOC_DF$GPS_Longitude,LOC_DF$GPS_Latitude), data=LOC_DF, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
kud_SPDF=spTransform(kud_SPDF,AKALBCRS)
kud_SPDF$AnNum=as.factor(kud_SPDF$AnNum)
kud_SPDF=kud_SPDF[,"AnNum"]

#kudmoose=kernelUD(kud_SPDF,kern="bivnorm",h="LSCV",hlim=c(0.001,1)) This does not work, as LSCV leads to very tiny values of smoothing (h). This is biologically unrealistic and is a severe downward bias on hr size.
kudmoose=kernelUD(kud_SPDF,kern="bivnorm",h="href") #Reference bandwidth tends to over-estimate HR size but the bias has emperically been found to be much smaller than lscv bias (See John Kie's 2013 paper). His ad hoc method is likely inappropriate due to the lacuna reducing rule (which in migratory moose would lead to extraordinary smoothing, and large hr estimates.)
kud_sizes=kernel.area(kudmoose, percent=seq(50, 99, by=1))

del_kud_size=kud_sizes
for (animal in 1:ncol(kud_sizes))
  {
  for (row in 1:nrow(kud_sizes))
  {
    ifelse((row==1),del_kud_size[row,animal]<-NA,del_kud_size[row,animal]<-(kud_sizes[row,animal]-kud_sizes[row-1,animal])/mean(kud_sizes[,animal]))
  }
  }
library(reshape)
df=melt(del_kud_size)
df$variable=rep(50:99,ncol(del_kud_size))
plot(df$variable,df$value) #All this work hammering away at this problem to basically go to 95%! But 95% is a data-justified cut-off in our case. ::shrug:: 

animal_kud_size=melt(kud_sizes["95",])
names(animal_kud_size)=c("AnNum","hectares")
animal_kud_size$AnNum=as.factor(substr(animal_kud_size$AnNum,2,8))

#Frustratingly, Plyr's ID function is masked by Adehabitat, so we need to unload and reload it.
detach(package:adehabitatHS, unload = TRUE)
detach(package:adehabitatHR, unload = TRUE)
detach(package:adehabitatLT, unload = TRUE)
detach(package:adehabitatMA, unload = TRUE)
detach(package:ade4, unload = TRUE)
detach(package:plyr, unload = TRUE)
library(plyr)

animal_kud_size$sex=join(animal_kud_size,masterid,by="AnNum",type="left",match="all")$Sex
animal_kud_size$CaptureGroup=substr(join(animal_kud_size,masterid,by="AnNum",type="left",match="all")$CapDate, nchar(as.character(join(animal_kud_size,masterid,by="AnNum",type="left",match="all")$CapDate))-4+1, nchar(as.character(join(animal_kud_size,masterid,by="AnNum",type="left",match="all")$CapDate)))#This code is literally satan. Did I really just use recursive joins?
#All this work and it turns out we only have 3 2015 animals with complete year of data. Therefore, the hr estimate is not biased low due to downriver HRs. Though, interestingly, 2 of the 3 2015 are positively tiny (12km2 and 40km2!)
boxplot((hectares*0.01)~sex,data=animal_kud_size,outline=FALSE)


#This code is because I frankly didn't believe the KD estimate of HR. They're so tiny! Much smaller than ballard found. Do we just have better tech than he did? Are we unlucky? Checking against MCP
library(ade4) #Reload all the parts I need.
library(adehabitatHR)
library(adehabitatHS)
library(adehabitatLT)
library(adehabitatMA)

mcpmoose=mcp(kud_SPDF, percent = 95)
#library(maptools)
#writePolyShape(mcpmoose, "homerange")
mcp_sizes=mcp.area(kud_SPDF, percent=seq(50, 100, by = 1))
del_mcp_size=mcp_sizes
for (animal in 1:ncol(mcp_sizes))
{
  for (row in 1:nrow(mcp_sizes))
  {
    ifelse((row==1),del_mcp_size[row,animal]<-NA,del_mcp_size[row,animal]<-(mcp_sizes[row,animal]-mcp_sizes[row-1,animal])/mean(mcp_sizes[,animal]))
  }
}

df2=melt(del_mcp_size)
df2$variable=rep(50:100,ncol(del_mcp_size))
plot(df2$variable,df2$value,ylim=c(0,2)) #Again, it turns out 94-95% is a pretty good break. At least this time I could c/p the code from before...

animal_mcp_size=melt(mcp_sizes["95",])
names(animal_mcp_size)=c("AnNum","hectares")
animal_mcp_size$AnNum=as.factor(substr(animal_mcp_size$AnNum,1,8))#For whatever reason, it didn't add an x to everything this time.

#Again with the loading and the unloading...
detach(package:adehabitatHS, unload = TRUE)
detach(package:adehabitatHR, unload = TRUE)
detach(package:adehabitatLT, unload = TRUE)
detach(package:adehabitatMA, unload = TRUE)
detach(package:ade4, unload = TRUE)
detach(package:plyr, unload = TRUE)
library(plyr)

animal_mcp_size$sex=join(animal_mcp_size,masterid,by="AnNum",type="left",match="all")$Sex
#At least I'm not doing the recrusive cohort matching code again.
par(mfrow=c(1,2))
boxplot((hectares*0.01)~sex,data=animal_mcp_size,outline=FALSE,main="MCP HR Est.", ylim=c(0,1200),ylab="HR area (square km)")
boxplot((hectares*0.01)~sex,data=animal_kud_size,outline=FALSE,main="KUD HR Est.", ylim=c(0,1200),ylab="HR area (square km)")
par(mfrow=c(1,1))
#Honestly, I'm pretty floored. MCP just about = KUD estimates. The moose really do have tiny home ranges.

quantile(subset(animal_kud_size,sex=="F")$hectare,c(0.25,0.5,0.75))*0.01
nrow(subset(animal_kud_size,sex=="F"))
quantile(subset(animal_kud_size,sex=="M")$hectare,c(0.25,0.5,0.75))*0.01
nrow(subset(animal_kud_size,sex=="M"))
boxplot((hectares*0.01)~sex,data=animal_kud_size,outline=FALSE, ylim=c(0,1200),ylab="Home range area (square km)", col=c("grey85","grey50"),names=c("Cow","Bull"))
mtext("Sex",side = 1,line = 2,cex = 1.25)
