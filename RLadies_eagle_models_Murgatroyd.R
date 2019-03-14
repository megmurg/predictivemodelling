
### R-Ladies Cape Town 14 March 2019 ####
# Predictive modelling
# Megan Murgatroyd | verreauxs@gmail.com | FitzPatrick Inst of African Ornitholgy

require(raster)
require(rgeos) #for gbuffer
require(spatialEco) #for vrm calculation
require(dismo) #for assessing model accuracy
#require(lme4) #for mixed models

#set working directory:
setwd("C:/Users/MRGSAR004/Documents/Verreauxs_eagles/Rladies")

#load GPS tracking data:
load("eaglesRLadies.RData")

#explore data
head(eagles)
str(eagles)

#plot data:
plot(eagles$longitude, eagles$latitude, col=eagles$name)

#set incorrect classes:
eagles$altitude=as.numeric(as.character(eagles$altitude))
eagles$speed_2d=as.numeric(as.character(eagles$speed_2d))
summary(eagles)

#add the nest location for each bird:
eagles$nest_lat=ifelse(eagles$name == "sto", -24.60695,
                       ifelse(eagles$name == "luc", -27.11595,
                              ifelse(eagles$name == "uil", -25.39721, "na")))
eagles$nest_lat=as.numeric(as.character(eagles$nest_lat))

eagles$nest_long=ifelse(eagles$name == "sto", 30.61312,
                        ifelse(eagles$name == "luc", 26.23938,
                               ifelse(eagles$name == "uil", 26.15159, "na")))
eagles$nest_long=as.numeric(as.character(eagles$nest_long))

#add nests to plot to check
points(eagles$nest_long, eagles$nest_lat, cex=0.8, pch=16, col="purple")

#### Pseudoabsence / Random points ####
#Create random points within a buffer for each animal
buff=12.73  #km

pseudo_e=NULL
IDs=unique(eagles$name)

for(i in 1:length(IDs)){
  bird=eagles[eagles$name == IDs[i],]             #subset data for one individual
  nest=bird[1:1,]                                 #extract first line of data for that id and call it nest
  coordinates(nest)=~nest_long + nest_lat         #nest the nest coordinate so it becomes a spatial point
  crs(nest)=crs("+proj=longlat +datum=WGS84 +zone=34s")   #set the current coordinate ref system
  nest=spTransform(nest, crs("+proj=utm +datum=WGS84 +zone=34s"))   #transform to UTM so we can measure buffer size in meters
  buffer=gBuffer(nest, width=buff*1000, byid=F)   #create a buffer (gBuffer is from rgeos package) in meters
  buffer=spTransform(buffer, crs("+proj=longlat +datum=WGS84 +zone=34s"))  #transform buffer back to longlat so the random points will be in the same coord-sys
  ps_points=spsample(buffer,n=(length(bird$name)*3),type="random")  #create random points. 3 times the number of real points per bird
  ps_points=as.data.frame(ps_points)              #transform to dataframe
  ps_points$name=c(as.character(bird$name))       #set class
  #head(ps_points)
  ps_points$nest_lat=bird$nest_lat[1]             #add nest lat and long
  ps_points$nest_long=bird$nest_long[1]
  names(ps_points)[names(ps_points) == "y"] <- "latitude"  #names to match real points dataframe
  names(ps_points)[names(ps_points) == "x"] <- "longitude"
  
  plot(ps_points$longitude, ps_points$latitude, col="red")
  points(bird$longitude, bird$latitude, col="blue")
  
  pseudo_e=rbind(pseudo_e, ps_points)
}

nrow(pseudo_e)/3==nrow(eagles) #TRUE  check there are the right number of points

#tidy:
rm(bird, buffer, nest, ps_points, buff, i, IDs)

#columns names need to match to be able to join dataframes:
names(pseudo_e)
names(eagles)

eagles$device_info_serial=NULL      #remove this column, it's not needed
pseudo_e$date_time="na"             #add date_time colum to random points and fill with "na"
eagles$date_time=as.character(eagles$date_time) #set as.character for the join
pseudo_e$altitude="na"
pseudo_e$speed_2d="na"

#add 'bin' column to both dataframes (will be used in analyses) and join
pseudo_e$bin=0
eagles$bin=1

eagles=rbind(eagles, pseudo_e)
summary(eagles)   #n=33272

#tidy:
rm(pseudo_e)

#### Annotate ####
# Attach terrain data to each location/observation

#load DEM for each individual
#SRTM 1arc-second (30m) available from https://earthexplorer.usgs.gov/
uil_dem<-raster("uil_dem.tif")
luc_dem<-raster("luc_dem.tif")
sto_dem<-raster("sto_dem.tif")

IDs=unique(eagles$name)
annotated=NULL

for(i in 1:length(IDs)){
  bird=eagles[eagles$name == IDs[i],]  #select bird
  coordinates(bird)=~longitude +latitude  #transform to spatial points
  
dem=if(bird$name[1]== "uil"){   #if statement to select correct DEM per bird
  uil_dem
} else if(bird$name[1]== "luc"){
  luc_dem
} else if(bird$name[1]== "sto"){
  sto_dem
} else {
  "na"
}

plot(dem)
points(bird)

slope=terrain(dem, opt=c('slope'), unit='degrees', neighbors=8)  #make slope layer from DEM
tri=terrain(dem, opt=c('tri'), neighbors=24)                     #make terrian ruggedness index
vrm=spatialEco::vrm(dem, s=5)                                    #make vector ruggedness measure

terrain.stack<-stack(list(elev=dem, slope=slope, tri=tri, vrm=vrm))  #create terrain stack for easy extractions

vars=extract(terrain.stack, bird)  #extract
vars=as.data.frame(vars)
bird=as.data.frame(bird)

bird=cbind(bird, vars)
annotated=rbind(annotated, bird)
}

#for the purpose of this example, remove NAs. 
#the nas for terrain data are caused by points off the raster extent
summary(annotated) #n=1
annotated<-annotated[complete.cases(annotated$vrm),]



#Add distance to nest using cartesian equation:
#This is needed for all points for the analyses
mdist=function(x,y,x2,y2){
  sqrt((x-x2)^2+(y-y2)^2)
}
#transorm to UTM to measure distance in meters:
eagle2=annotated #make a copy of the latlong data
coordinates(eagle2)=~longitude + latitude  #project fix locations
crs(eagle2)=CRS("+proj=longlat +datum=WGS84 +zone=34s") #give current projection
eagle2=spTransform(eagle2, CRS("+proj=utm +datum=WGS84 +zone=34S")) #transform to UTM
eagle2=as.data.frame(eagle2) #return to dataframe
#same as above for nest location:
coordinates(eagle2)=~nest_long + nest_lat  #transform nest locations
crs(eagle2)=CRS("+proj=longlat +datum=WGS84 +zone=34s")
eagle2=spTransform(eagle2, CRS("+proj=utm +datum=WGS84 +zone=34S"))
eagle2=as.data.frame(eagle2)
head(eagle2) #fix and nest location in utm

#apply mdist:
annotated$nest_dist=mdist(eagle2$longitude, eagle2$latitude, eagle2$nest_long, eagle2$nest_lat)/1000
summary(annotated$nest_dist)

# #another function to measure distance:
# require(geosphere)
# annotated$nest_dist=distHaversine(annotated[,c("latitude","longitude")], annotated[,c("nest_lat","nest_long")], r=6378137)/1000
# summary(annotated$nest_dist)
# summary(annotated)

#### Classfiying RISK ####
# Wind turbine collision risk is any point < 200m above ground level
head(annotated)
annotated$altitude=as.numeric(as.character(annotated$altitude)) #set class - this is altiude above sea level
annotated$agl=annotated$altitude-annotated$elev #calculate altitude above ground level
annotated$agl=as.numeric(as.character(annotated$agl)) #set class
summary(annotated$agl) #NAs are the pseudo points sinces they don't have an initial altitude

#remove non-risky points and match length 3xpseudo points:
bin1=subset(annotated, bin==1)
bin0=subset(annotated, bin==0)

#remove negative agl values:
bin1=bin1[(bin1$agl >= 0),]   #n 8317 > 6536

#remove non-risky values (above 200m):
bin1=bin1[(bin1$agl <= 200),]   #n 6536 > 5493

#match 3xrisk points and rejoin
#should really be per individual... but anyway...
bin0=bin0[sample(nrow(bin0), nrow(bin1)*3), ]
annotated=rbind(bin1, bin0)
rm(bin1, bin0)

#tidy:
rm(bird, dem, eagle2, luc_dem, terrain.stack, tri, uil_dem, vars, vrm, i, IDs, mdist, slope, sto_dem, eagles)

### ANALYSES ####
#spilt data to train and test set.
test=subset(annotated, name=="luc")
train=subset(annotated, name!="luc")

#check for correlations # >0.6 = too correlated for same model.
str(annotated)
check=annotated

check$name=NULL
check$date_time=NULL
check$speed_2d=NULL

cor(check)
rm(check)

#simple model:
model=glm(bin~ scale(nest_dist) + scale(elev)+ scale(vrm) + scale(slope), data=train, family=binomial, na.action = na.fail)
summary(model)

# #with random effect
# require(lme4)
# model=glmer(bin~ scale(nest_dist) + scale(elev)+ scale(vrm) + scale(slope)+(1|name), data=train, family=binomial, na.action = na.fail)
# summary(model)


### PREDICTION ####
#explore data
head(test)
str(test)

#load raster
dem=raster("luc_dem.tif")
plot(dem)

#make rasterbrick
slope=terrain(dem, opt=c('slope'), unit='degrees', neighbors=8)
vrm=spatialEco::vrm(dem, s=5) #s=odd number for grid calculation e.g. 3x3 or 5x5...

terrain.stack<-stack(list(elev=dem, slope=slope, vrm=vrm))
plot(terrain.stack)

#Crop terrain stack by buffer around nest:
buff=12.73

nest=cbind(nest_lat=c(-27.11595), nest_long=c(26.23938)) #nest location
nest=as.data.frame(nest)
coordinates(nest)=~nest_long + nest_lat
crs(nest)=crs("+proj=longlat +datum=WGS84 +zone=34s") #project
nest=spTransform(nest, crs("+proj=utm +datum=WGS84 +zone=34s")) #transform to UTM
buffer=gBuffer(nest, width=buff*1000, byid=F) #make buffer
buffer=spTransform(buffer, crs("+proj=longlat +datum=WGS84 +zone=34s")) #transform back
terrain.stack=crop(terrain.stack, buffer) #crop terrain stack by buffer extents
#terrain.stack=mask(terrain.stack, buffer) #use this if you want to make it circular

rb_pts=rasterToPoints(terrain.stack)  #transform terrain.stack to dataframes with col for each variable
rb_pts=as.data.frame(rb_pts)
head(rb_pts)

# x=rb_pts$x
# y=rb_pts$y

rb_pts$nest_lat= -27.11595 
rb_pts$nest_long= 26.23938

#tranform to utm and measure distances to nest:
coordinates(rb_pts)=~nest_long + nest_lat  #project & transform fix locations
# plot(dem)
# points(rb_pts[1,])
crs(rb_pts)=CRS("+proj=longlat +datum=WGS84 +zone=34S")
rb_pts=spTransform(rb_pts, CRS("+proj=utm +datum=WGS84 +zone=34S"))
rb_pts=as.data.frame(rb_pts)

coordinates(rb_pts)=~x + y  #transform nest locations
crs(rb_pts)=CRS("+proj=longlat +datum=WGS84")
rb_pts=spTransform(rb_pts, CRS("+proj=utm +datum=WGS84 +zone=34S"))
rb_pts=as.data.frame(rb_pts)

#measure nest_dist:
mdist=function(x,y,x2,y2){
  sqrt((x-x2)^2+(y-y2)^2)
}

rb_pts$nest_dist=mdist(rb_pts$nest_long, rb_pts$nest_lat, rb_pts$x, rb_pts$y)/1000

#transform x y back to longlat:
coordinates(rb_pts)=~x + y  #transform nest locations
crs(rb_pts)=CRS("+proj=utm +datum=WGS84 +zone=34S")
rb_pts=spTransform(rb_pts, CRS("+proj=longlat +datum=WGS84"))
rb_pts=as.data.frame(rb_pts)
summary(rb_pts)

#PREDICTED RISK PLOT:
pred=predict(model, rb_pts, re.form = NA, type = "response", na.action = na.fail)
pred=as.data.frame(pred)
summary(pred$pred)
head(rb_pts)

#bind the data needed to make a risk plot
toplot=cbind(x=c(rb_pts$x), y=c(rb_pts$y), pred=c(pred$pred))
toplot=as.data.frame(toplot)
head(toplot)

# create spatial points data frame
coordinates(toplot) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(toplot) <- TRUE
# coerce to raster
rasterDF <- raster(toplot)
rasterDF
plot(rasterDF)
plot(buffer, add=T)

#add test "known" risk data to plot:
#only need bin ==1
test1=subset(test, bin==1)
coordinates(test1) <- ~ longitude + latitude
points(test1, cex=0.5)

#### Assess model accuracy ####
coordinates(test) <- ~ longitude + latitude
test_pred=extract(rasterDF, test)
test=as.data.frame(test)
test$pred=test_pred

require(dismo)
e = evaluate(p=subset(test, bin==1)$pred, a=subset(test, bin==0)$pred)
e

