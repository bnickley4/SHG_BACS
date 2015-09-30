###BACS 9-26-15 BN###

###read in BACS data and burn data
dat=read.csv("Rdat_BN.csv")
burn=read.csv("burn.vectors.csv")

###subset data 
j<-dat[which(dat$PointName=='J1'),]	#BACS data (dat) by point J1
head(j)
j.day=j$DaysAfterN
j.burn=burn[,1]	#burn dates for J1
j.burn<-j.burn[!is.na(j.burn)]	#remove na's

vars=c("PointName","DaysAfterN")	#make a vector w/ variable names
dat.ss<-dat[vars]	#subset data including only variables of interest
head(dat.ss)

###get days-since-burned (dsb)
j.dsb<-NA	#create an empty object
for(i in 1:length(j.day)){	#for each value in j.day
for(k in 1:length(j.burn)){	#and for each value in j.burn
if(j.burn[k]<=j.day[i]){	#if a value in j.burn is less than any j.day value
j.dsb[i]<-j.day[i]-j.burn[k]}}}	#then, take the value of j.day and subtract it from the value of j.burn
	#and store the values in a vector with length = length(j.day)
	#R goes through the code sequentially, so b/c j.burn days in order (sm. to lg.), only the max value of j.burn that is <= j.day is stored in j.delta, thankfully
j.dsb

###obtain burn date and survey date vectors for each point
p.name<-levels(dat$PointName)	#point names in vector "p.name"
dat.s<-split(dat.s, dat.s$PointName)	#split by pt. name, **see 'subset data' above
head(dat.s$J1[,2])	#vector of dates for J1
#library(plyr)	#recode variable 'PointName' into integers
#dat.ss$pcode<-mapvalues(dat.ss$PointName, from=p.name, to=c(1:length(p.name)))
#head(dat.ss)
#vars=c("pcode", "DaysAfterN")
#dat.ss<-dat.ss[vars]
#head(dat.ss)
#dat.s<-split(dat.ss, dat.ss$PointName)
#colnames(burn)<-c(1:length(colnames(burn)))	#Don't need to do this#recode point names of burn into integers as well

###########not necessary##########
dat.s$J16[,2]	###=dat.s[[p]][,2] in loop###vector of survey dates for point 73
burn[,3]	#vector of burn dates for point 73
p.name[3]	#"J16" is point 3

#for(p in 1:ncol(burn)){for(i in 1:length(dat.s[[p]][,2])){for(k in 1:length(burn[,p])){if(burn[k,p]<=dat.s[[p]][i,2]){dsb[i]<-dat.s[[p]][i,2]-burn[k,p]}}}}
##getting an error here, prob. b/c of NAs in 'burn'
#max(dat$DaysAfterN) #=5972, this is the last survey date. i can set burn NAs to >6000 and they will be ignored
burn[is.na(burn)] <- 9999	#replace NAs w/ 9999

###blank matrix
dsb<-matrix(data="NA", ncol=79, nrow=50)	#dsb "days since burned"

###The Magic Code refined- get days since burned!###
for(p in 1:ncol(burn)){for(i in 1:length(dat.s[[p]][,2])){for(k in 1:length(burn[,p])){if(burn[k,p]<=dat.s[[p]][i,2]){dsb[i,p]<-dat.s[[p]][i,2]-burn[k,p]}}}}
#same code broken down
for(p in 1:ncol(burn)){	#burn has 79 columns, 1:79 corresponding to points J1:R8
for(i in 1:length(dat.s[[p]][,2])){	#dat.s is my survey dates split by point (1:79). I'm selecting split/point 'p' column 2 [,2] has the dates
for(k in 1:length(burn[,p])){	#selecting column 'p' from burn which cooresponds to point (1:79)
if(burn[k,p]<=dat.s[[p]][i,2]){	#if the value(date) 'k' from the burn vector 'p' <= value (date) 'i' from dat.s vector 
dsb[i,p]<-dat.s[[p]][i,2]-burn[k,p]}}}}	#write the difference of the two into my blank matrix 'dsb'

###validation
p.name[73]	#PointName 73 is pcode L5
day.v<-dat.s$Q2[,2]	#make a vector for survey dates at point '73' = 'L5'
burn.v<-burn[,73]	#make a vector for burn dates for point '73'
dsb.v<-NA	#create blank object 'dsb.v'
for(i in 1:length(day.v)){for(k in 1:length(burn.v)){if(burn.v[k]<=day.v[i]){dsb.v[i]<-day.v[i]-burn.v[k]}}}
#(above)use old code to create days-since-burned vector just for point 73
dsb.v
# [1] 399 401 409  35  37  43 406 408 414   3   9  14

dsb[,73]
# [1] "399" "401" "409" "35"  "37"  "43"  "406" "408" "414" "3"   "9"   "14"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA" 
#[20] "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA" 
#[39] "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA"  "NA" 
#cooresponds with days-since-burned vector generated w/ code ran to create 'dsb' matrix containig days-since-burned columns (1:79) for each point.
colnames(dsb)<-names(burn)
#write.csv(dsb, file = "dsb.csv")	#write csv of 'dsb' data frame to working directory (R data)