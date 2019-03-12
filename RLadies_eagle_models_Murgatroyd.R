
### R-Ladies ####
# Predictive modelling : March 2019 
# Megan Murgatroyd | verreauxs@gmail.com | FitzPatrick Inst of African Ornitholgy

#set working directory:
setwd("C:/Users/MRGSAR004/Documents/Verreauxs_eagles/Rladies")
                  
#load tracking data:
load("eagle2_RLadies_Murgatroyd.RData")

#explore data
head(eagle2)
str(eagle2)

#set incorrect classes:
eagle2$altitude=as.numeric(as.character(eagle2$altitude))
eagle2$speed_2d=as.numeric(as.character(eagle2$speed_2d))
summary(eagle2)

#### Pseudoabsence / Random points ####
buff=12.73277

#nest location:
nest=NULL
nest$nest_lat= -18.60695
nest$nest_long= 36.61312
nest=as.data.frame(nest)
coordinates(nest)=~nest_long +nest_lat

#transform to UTM
crs(nest)=CRS("+proj=longlat +datum=WGS84")
nest=spTransform(nest, CRS("+proj=utm +datum=WGS84 +zone=34S"))

#12km nest buffer:
buffer=gBuffer(nest, width=(buff*1000), byid=F)
buffer=spTransform(buffer, CRS("+proj=longlat +datum=WGS84"))

#generate 3 times random points inside buffer:
set.seed(10)
ps_points=spsample(buffer,n=(length(eagle2$name)*3),type="random")

#check:
plot(buffer)
points(ps_points, cex=0.3, col="blue")

#make columns match and join dataframes:
ps_points=as.data.frame(ps_points)
ps_points$name=c(as.character(eagle2$name))

ps_points$name=c(as.character(eagle2$name))
ps_points$device_info_serial=c(as.character(eagle2$device_info_serial))
ps_points$date_time="na"
eagle2$date_time=as.character(eagle2$date_time)
ps_points$altitude="na"
ps_points$speed_2d="na"
colnames(ps_points)[colnames(ps_points)=="x"] <- "longitude"
colnames(ps_points)[colnames(ps_points)=="y"] <- "latitude"

#add bin and join
ps_points$bin=0
eagle2$bin=1

eagle2=rbind(eagle2, ps_points)

#tidy:
rm(buff, buffer, nest, ps_points)

#### Annotate ####
# Attach terrain data to each location/observation
require(raster)
#load rasterBrick
rb2=stack("rbrick2.tif")

#name brick layers:
names(rb2)[1] <- 'aspect'
names(rb2)[2] <- 'slope'
names(rb2)[3] <- 'slope_sd3'
names(rb2)[4] <- 'slope_sd5'
names(rb2)[5] <- 'alt_sd3'
names(rb2)[6] <- 'alt_sd5'
names(rb2)[7] <- 'alt'
names(rb2)[8] <- 'tri'
names(rb2)[9] <- 'vrm'

plot(rb2)

#transform eagle data to spatial points:
coordinates(eagle2) <- ~ longitude + latitude
crs(eagle2)=CRS("+proj=longlat +datum=WGS84")

#visualise:
plot(rb2$alt)
points(eagle2, cex=0.5, col=eagle2$bin)

#add nest:
nest=NULL
nest$nest_lat= -18.60695
nest$nest_long= 36.61312
nest=as.data.frame(nest)
coordinates(nest)=~nest_long +nest_lat
points(nest, pch=16, col="red")

#extract data:
vars<-extract(rb2, eagle2)
eagle2=as.data.frame(eagle2)
vars=as.data.frame(vars)
eagle2=cbind(eagle2, vars)
summary(eagle2)

#remove spatial outliers
eagle2=subset(eagle2, complete.cases(vrm))

#add distance to nest for each fix:
nest=as.data.frame(nest)
head(nest)
eagle2$nest_lat=nest$nest_lat
eagle2$nest_long=nest$nest_long

#function for calc distance in m (from utm):
mdist=function(x,y,x2,y2){
  sqrt((x-x2)^2+(y-y2)^2)
}
#tranform to utm and measure distances to nest:
eagle2_ll=eagle2 #male a copy of the latlong data
coordinates(eagle2)=~longitude + latitude  #project & transform fix locations
crs(eagle2)=CRS("+proj=longlat +datum=WGS84") #give current projection
eagle2=spTransform(eagle2, CRS("+proj=utm +datum=WGS84 +zone=34S")) #transform to UTM
eagle2=as.data.frame(eagle2) #return to dataframe

#same as above for nest location:
coordinates(eagle2)=~nest_long + nest_lat  #transform nest locations
crs(eagle2)=CRS("+proj=longlat +datum=WGS84")
eagle2=spTransform(eagle2, CRS("+proj=utm +datum=WGS84 +zone=34S"))
eagle2=as.data.frame(eagle2)
head(eagle2) #fix and nest loca in utm

#apply mdist:
eagle2_ll$nest_dist=mdist(eagle2$longitude, eagle2$latitude, eagle2$nest_long, eagle2$nest_lat)/1000
summary(eagle2_ll$nest_dist)

#tidy:
eagle2=eagle2_ll
rm(eagle2_ll, nest, vars, mdist, rb2)

#### Classfiying RISK ####
# Wind turbine collision risk is any point < 200m above ground level
eagle2$altitude=as.numeric(as.character(eagle2$altitude))
eagle2$agl=eagle2$altitude-eagle2$alt
eagle2$agl=as.numeric(as.character(eagle2$agl))
summary(eagle2$agl)

#remove non-risky points and match length 3xpseudo points:
bin1=subset(eagle2, bin==1)
bin0=subset(eagle2, bin==0)

#remove negative values:
bin1=bin1[(bin1$agl >= 0),]

bin1$risk=ifelse(bin1$agl <= 230, 1, 0)
bin1$risk=as.factor(bin1$risk)
summary(bin1$risk)

#remove non-risky values:
bin1=bin1[(bin1$risk == 1),]
bin1$risk=NULL

#match 3xrisk points and rejoin
bin0=bin0[sample(nrow(bin0), nrow(bin1)*3), ]
eagle2=rbind(bin1, bin0)
rm(bin1, bin0)

### ANALYSES ####
#check for correlations # >0.6 = too correlated for same model.
str(eagle2)
check=eagle2

check$name=NULL
check$device_info_serial=NULL
check$date_time=NULL
check$altitude=NULL
check$speed_2d=NULL

cor(check)
rm(check)

#simple model:
model=glm(bin~ scale(nest_dist) + scale(alt)+ scale(vrm) + scale(slope), data=eagle2, family=binomial, na.action = na.fail)
summary(model)


### PREDICTION ####
#load tracking data:
load("eagle1_RLadies_Murgatroyd.RData")

#explore data
head(eagle1)
str(eagle1)

#set incorrect classes:
eagle1$altitude=as.numeric(as.character(eagle1$altitude))
eagle1$speed_2d=as.numeric(as.character(eagle1$speed_2d))
summary(eagle1)

#load rasterBrick
rb1=stack("rbrick1.tif")

#name brick layers:
names(rb1)[1] <- 'aspect'
names(rb1)[2] <- 'slope'
names(rb1)[3] <- 'slope_sd3'
names(rb1)[4] <- 'slope_sd5'
names(rb1)[5] <- 'alt_sd3'
names(rb1)[6] <- 'alt_sd5'
names(rb1)[7] <- 'alt'
names(rb1)[8] <- 'tri'
names(rb1)[9] <- 'vrm'

rb_pts=rasterToPoints(rb1)
rb_pts=as.data.frame(rb_pts)
head(rb_pts)

x=rb_pts$x
y=rb_pts$y

rb_pts$nest_lat= -19.39721 
rb_pts$nest_long= 32.15159

#tranform to utm and measure distances to nest:
coordinates(rb_pts)=~nest_long + nest_lat  #project & transform fix locations
crs(rb_pts)=CRS("+proj=longlat +datum=WGS84")
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

#transform back x y:
coordinates(rb_pts)=~x + y  #transform nest locations
crs(rb_pts)=CRS("+proj=utm +datum=WGS84 +zone=34S")
rb_pts=spTransform(rb_pts, CRS("+proj=longlat +datum=WGS84"))
rb_pts=as.data.frame(rb_pts)
summary(rb_pts)

#PREDICTED RISK PLOT:
pred=predict(model, rb_pts, re.form = NA, type = "response", na.action = na.fail)
pred=as.data.frame(pred)
summary(pred$pred)

toplot=cbind(x=c(rb_pts$x), y=c(rb_pts$y), pred=c(pred$pred))
head(toplot)
toplot=as.data.frame(toplot)
head(toplot)

# create spatial points data frame
coordinates(toplot) <- ~ x + y
#toplot=rasterize(toplot)
# coerce to SpatialPixelsDataFrame
gridded(toplot) <- TRUE
# coerce to raster
rasterDF <- raster(toplot)
rasterDF
plot(rasterDF)

####eagle risk points:####

coordinates(eagle1)=~longitude + latitude  #project & transform fix locations
crs(eagle1)=CRS("+proj=longlat +datum=WGS84")

#extract data:
vars<-extract(rb1, eagle1)
eagle1=as.data.frame(eagle1)
vars=as.data.frame(vars)
eagle1=cbind(eagle1, vars)
summary(eagle1)

# Wind turbine collision risk is any point < 200m above ground level
eagle1$altitude=as.numeric(as.character(eagle1$altitude))
eagle1$agl=eagle1$altitude-eagle1$alt
eagle1$agl=as.numeric(as.character(eagle1$agl))
summary(eagle1$agl)

eagle1$risk=ifelse(eagle1$agl <= 230, 1, 0)
eagle1$risk=as.factor(eagle1$risk)
summary(eagle1$risk)

#remove non-risky values:
eagle1=eagle1[(eagle1$risk == 1),]
eagle1=eagle1[complete.cases(eagle1$latitude),]
eagle1$risk=NULL

summary(eagle1)
coordinates(eagle1)=~longitude + latitude  #project & transform fix locations
crs(eagle1)=CRS("+proj=longlat +datum=WGS84")
points(eagle1, cex=0.2)


