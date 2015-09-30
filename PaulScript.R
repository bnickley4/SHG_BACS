#The following script file tracks changes make to master 
#Sandhills Game Lands (SHGL)BACS data to format for analysis in
#UNMARKED.  

#     Created: 5/26/14
#        By Paul Taillie

# Set working directory (change if working from a computer other 
#   than  WRC #26220
setwd("C:/Users/tailliep/Documents/BACS/SHGL_longterm")

# Load required packages
library(unmarked)
#############################################################
#  subset data

# Load data
dat<-read.csv("SHGL_Unmarked2014.csv")
head(dat)

# Get list of points that have been surveyed at least 20 times
freq<-as.data.frame(table(dat$PointName)) #creates table of number of surveys
freq<-freq[order(freq$Freq),]#sort by number of times surveyed
freq.sub<-subset(freq, freq$Freq>20)
points<-as.vector(freq.sub[,1])

#subset original data
keep<-dat$PointName %in% points
dat<-dat[keep==TRUE,]

dat<-dat[dat$SurveyID!=4,]
dat<-dat[dat$SurveyID!=5,]

#############################################################
#    Convert to wide format

dat.in<-data.frame(dat$PointName,as.integer(dat$SurveyID),dat$TotalBACS,dat$Temp,dat$Wind,dat$Cloud,dat$Date)
colnames(dat.in)<-c("Point","Survey","Count","Temp","Wind","Cloud","Date")

# Fix wind values
a<-which(dat.in$Wind=="0-10")
dat.in$Wind[a]<-5
a<-which(dat.in$Wind=="0-")
dat.in$Wind[a]<-0
a<-which(dat.in$Wind=="0-20")
dat.in$Wind[a]<-10
a<-which(dat.in$Wind=="0-15")
dat.in$Wind[a]<-5
a<-which(dat.in$Wind=="0-3")
dat.in$Wind[a]<-2
a<-which(dat.in$Wind=="0-5")
dat.in$Wind[a]<-3
a<-which(dat.in$Wind=="15-May")
dat.in$Wind[a]<-10
a<-which(dat.in$Wind=="20-May")
dat.in$Wind[a]<-15
dat.in$Wind<-as.integer(dat.in$Wind)

# remove junk
a<-which(dat.in$Count!="NA")
dat.in<-dat.in[a,]

#Put in wide format
# first make some empty data frames
ydat<-matrix(data = NA, nrow = 80, ncol = 29)
ydat<-as.data.frame(ydat)
rownames(ydat)<-unique(dat.in$Point)
temp<-matrix(data = NA, nrow = 80, ncol = 29)
temp<-as.data.frame(temp)
rownames(temp)<-unique(dat.in$Point)
wind<-matrix(data = NA, nrow = 80, ncol = 29)
wind<-as.data.frame(wind)
rownames(wind)<-unique(dat.in$Point)
cloud<-matrix(data = NA, nrow = 80, ncol = 29)
cloud<-as.data.frame(cloud)
rownames(cloud)<-unique(dat.in$Point)

for (i in 1:nrow(dat.in)){
tempsite<-as.character(dat.in$Point[i])
ydat[tempsite,dat.in$Survey[i]]<-dat.in$Count[i]
temp[tempsite,dat.in$Survey[i]]<-dat.in$Temp[i]
wind[tempsite,dat.in$Survey[i]]<-dat.in$Wind[i]
cloud[tempsite,dat.in$Survey[i]]<-dat.in$Cloud[i]
}

#remove SurveyID 4 and 5 (to keep the secondary sampling periods
#  to three for every year
ydat$V4<-NULL
ydat$V5<-NULL
temp$V4<-NULL
temp$V5<-NULL
wind$V4<-NULL
wind$V5<-NULL
cloud$V4<-NULL
cloud$V5<-NULL

# standardize covariates
temp=as.matrix(temp)
wind=as.matrix(wind)
cloud=as.matrix(cloud)
temp.mean=mean(c(temp),na.rm=T)
temp.sd=sd(c(temp),na.rm=T)
temp=(temp-temp.mean)/temp.sd
wind.mean=mean(c(wind),na.rm=T)
wind.sd=sd(c(wind),na.rm=T)
wind=(wind-wind.mean)/wind.sd
cloud.mean=mean(c(cloud),na.rm=T)
cloud.sd=sd(c(cloud),na.rm=T)
cloud=(cloud-cloud.mean)/cloud.sd


umf<-unmarkedFramePCO(y=ydat,
                      obsCovs=list(temp=as.matrix(temp),
					     wind=as.matrix(wind),
					     cloud=as.matrix(cloud)),
				numPrimary=9)
mod1<-pcountOpen(~1,   ~1,  ~1,   ~temp+wind+cloud,data=umf)
mod2<-pcountOpen(~1,   ~1,  ~1,   ~wind+cloud,data=umf)


# Reverse transform parameters
lam <- exp(coef(mod1, type="lambda"))
gam <- exp(coef(mod1, type="gamma"))
om <- plogis(coef(mod1, type="omega"))
p <- plogis(coef(mod1, type="det"))


# write a new csv of formatted table for use in UNMARKED
write.csv(



