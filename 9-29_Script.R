###BACS 9-28-15 BN###

#Trouble shooting 'dsb' data#
#dsb = days since burn

dat=read.csv("Rdat_BN.csv") #load table of all data for all years
burn=read.csv("burn.vectors.csv") #each column is a point, rows for day burned (in days since 1999)

unique(dat$PointName)
p.name<-levels(dat$PointName)
p.names<-p.name
p.vars<-names(burn)

##these should both have 79 unique point names.  dat has data for many points not run in 2015
##the code below finds point names that differ from the 2015 point names (names in burn)
extra.points<-vector(mode="numeric",length=200)	#make empty vector
for(i in 1:length(p.vars)){for(k in 1:length(p.names)){if(p.vars[i]==p.names[k]){p.names[k]<-"NA"}}}	#get point names that do NOT overlap p.vars (the 79 point names on the maps, run in 2015)
p.names					####point names that do NOT overlaps w/ 2015 point names or burn vectors####
p.names<-p.names[p.names!="NA"]	#remove "NA"'s
write.csv(p.names, "p.names.csv")	#write a csv w/ this data
length(p.name)-length(p.names)	#this equals 79 as expected

dat.r<-dat	#copy of dat

point.v<-vector(mode="character", length=length(dat.r$PointName))	#make a blank vector of length dat.r$PointName, mode="character" important, modes must match!!!!
point.v<-dat.r$PointName[1:length(dat.r$PointName)]	#put data from data frame dat.r into vector for PointName
head(point.v)
point.vv<-point.v	#a copy of point.v that will be used when replacing values in the loop

#make a loop to replace non-2015 point names w/ NAs in the vector 'point.vv'
for(i in 1:length(point.v)){for(k in 1:length(p.names)){if(p.names[k]==point.v[i]){point.vv[i]<-"NA"}}}
dat.r$PointName<-as.character(dat.r$PointName)				#make point name a character
dat.r$PointName<-point.vv							#replace value in dat.r$PointName w/ values from vector point.vv generated in loop (above)
unique(dat.r$PointName)								#look at unique values. notice, there is an "NA" in there which is what replaced point names that do not coorespond w/ a burn vector (i.e. not run in 2015, not on SHG maps BN made)
length(unique(dat.r$PointName))						#length is 80, 79 + the "NA"
length(which(dat.r$PointName=="NA"))					#how many "NA"'s did I generate?
length(dat.r$PointName)								#how many total observations
length(dat.r$PointName)-length(which(dat.r$PointName=="NA"))	#difference=num. of non-"NA"'s left

#turn "NA"'s into real NAs...
length(dat.r$PointName[!is.na(dat.r$PointName)])	#the "NA"'s are NOT being considered as NAs by R
dat.r$PointName[dat.r$PointName=="NA"]=NA			#make them real NAs
length(dat.r$PointName[!is.na(dat.r$PointName)])	#it worked!

#remove rows w/ NA's
dat.omit<-na.omit(dat.r)	#make an new data frame w/ NAs excluded
length(dat.omit$PointName)	#check the length, it checks out.

###okay, but i really need to find out what is going on w/ the names and why names from burn vectors, maps and 2015 survey do not coorespond 100% w/ names from previous years.
###i'm not just going to delete half of our obs. i bet those 130 points w/ diff. names are 'misnamed' but actually coorespond to one of the 1:79 point names

#subset data
names(dat.omit)
vars<-names(dat.omit)[c(6,4,19,9,11,10,12:15)]	#define variables to use in unmarked
vars
dat.t<-dat.omit[vars]		#subset dat.omit keeping the variables defined in 'vars'


####then i made the burn vectors following the script from 9-26####

dsb.df<-read.csv("dsb_NEWu.csv")	#read in the newest dsb vectors

#now I want to take the burn vectors and put them into the dat.n data frame...
dsb<-c(rep(NA,length(dat.t$PointName)))
dsb<-as.vector(dsb, mode="integer")	#make an empty vector to add to dat.t


dat.t$new.col<-dsb				#add the vector to dat.t
which(names(dat.t)=="new.col")		#get position of the new column
names(dat.t)[12]<-"dsb"				#rename new.col to dsb

dat.n<-split(dat.t, dat.t$PointName)	#split data by 'pcode' which is the point #'s
dat.nc<-dat.n				#make a copy of dat.n for later use (not really necessary)

dat.n$J12[1,12]	#here's the first value for the new column 'dsb'
dsb.df$J12[1]	#and here's the value I want to put in there...

###this code puts the 'dsb' data into the dat.nc dataframe###
for(i in 1:79){for(k in 1:length(dat.n[[i]][,12])){if(dsb.df[[i]][k]!="NA"){dat.nc[[i]][k,12]<-dsb.df[[i]][k]}}}
dat.nc$J12[1,12]	#value placed in here is from 'dsb'
dat.n$J12[1,12]	#still an "NA" in dat.n

###unsplit 
dat.nc<-do.call("rbind", dat.nc)	#I don't know why 'unsplit' doesn't work...

length(which(is.na(dat.nc$dsb)))	#any NAs left? NO
#write.csv(dat.nc, file="BACS_dsb.csv")	#write the .csv

########################################################################################
###"UNK"'s did not get coded into the dataframe :(
dat.nc$dsb[c(which(dat.nc$pcode=="Q2"))]
#it is off by 3 places (61 should be 59), or something.  indexing wrong.

#ANSWER -> > dat.nc[64] 	#ORDER NUMBER BEING INDEXED, NOT PCODE NUMBER
#$`67`
#head(dat.nc$Q2)		#order shifted when split, 1,10,11... not 1,2,3...
#this means the dsb vectors are WRONG TOO for the same reason.
#don't use pcode. only order matters. split by PointName

#clean-up
rm(dat.n)
rm(dat.nc)
rm(dsb)
rm(dsb.df)
