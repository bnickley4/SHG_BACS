###################################################################
### Code for module 8: Introduction to the modeling of occurrence 
### and species distributions
###################################################################


### 8.4. Simulation of data and first analysis using unmarked and WinBUGS/JAGS

# Create a covariate called vegHt
nSites <- 100
set.seed(443)                   # so that we all get the same values of vegHt
vegHt <- sort(runif(nSites, 1, 3)) # sort for convenience

# Suppose that occupancy probability increases with vegHt
# The relationship is described by an intercept of -3 and
#    a slope parameter of 2 on the logit scale
psi <- plogis(-3 + 2*vegHt)

# Now we go to 100 sites and observe their occurrence strate (perfectly)
z <- rbinom(nSites, 1, psi)

# We can fit a model that relates abundance to vegHt using the glm() function
#  with "family=binomial"
summary(fm.glm1 <- glm(z ~ vegHt, family=binomial))

# Do some analysis of the results
plot(vegHt, z, xlab="Vegetation height", ylab="Occurrence (z)")
glm1.est <- coef(fm.glm1)
plot(function(x) plogis(-3 + 2*x), 1, 3, add=TRUE, lwd=3)
plot(function(x) plogis(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,
     lwd=3, col="blue")
legend(2.5, 0.2, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=3)

nVisits <- 3
wind <- array(rnorm(nSites * nVisits), dim = c(nSites, nVisits))
p <- plogis(1 - 2*wind)
# plot(p ~ wind)
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(nVisits, z[i], p[i,])
}

# Look at the data
cbind(z=z, y1=y[,1], y2=y[,2], y3=y[,3])

# Load library, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(vegHt = vegHt), obsCovs = list(wind = wind))
summary(umf)

# Fit a model and extract estimates
# Detection covariates follow first tilde, then come occupancy covariates
# We don’t time run: occupancy models much faster than Nmix models !
summary(fm.occ <- occu(~wind ~vegHt, data=umf))

# All estimates are on logit-scale
par(mfrow = c(1,2))
beta1 <- coef(fm.occ) 
plot(function(x) plogis(beta1[1] + beta1[2]*x), 1, 3, xlab="Vegetation height", ylab="Occupancy probability", ylim = c(0, 1))
plot(function(x) plogis(beta1[3] + beta1[4]*x), -3, 3, xlab="Wind", ylab="Detection probability", ylim = c(0, 1))

# Or predictions for new values of vegHt, say 1.2 and 3.1
newdat <- data.frame(vegHt=c(1.2, 3.1))
predict(fm.occ, type="state", newdata=newdat)

ranef(fm.occ)

sum(ranef(fm.occ)@post[,2,])

# Bundle data
win.data <- list(y = y, vegHt = vegHt, wind = wind, R = nrow(y), T = ncol(y))

# Specify model in BUGS language
sink("model.txt")
cat("
model {

# Priors
alpha.occ ~ dunif(-10, 10)
beta.occ ~ dunif(-10, 10)
alpha.p ~ dunif(-10, 10)
beta.p ~ dunif(-10, 10)

# Likelihood
for (i in 1:R) {
   # True state model for the partially observed true state
   z[i] ~ dbern(psi[i])             # True occupancy z at site i
   logit(psi[i]) <- alpha.occ + beta.occ * vegHt[i]

   for (j in 1:T) {
      # Observation model for the actual observations
      y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p + beta.p * wind[i,j]
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])       # Number of occupied sites among those studied
}
",fill = TRUE)
sink()


# Initial values
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.occ", "beta.occ", "alpha.p", "beta.p", "occ.fs")
# add "z" for estimates of the random effects

# MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 3

# Call WinBUGS from R (ART 0.28 min)
out1 <- bugs(win.data, inits, params, "model.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Call JAGS from R (ART 0.26 min)
library("R2jags")		# requires rjags
system.time(out2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb) )
traceplot(out2)

# Summarize posteriors
print(out2, dig = 3)

# Graphical comparison between truth and estimates
par(mfrow = c(1,2), mar = c(5,4,4,4), cex.lab = 1.2, cex.main = 1.2)
predHt <- seq(1,3,,100)
predW <- seq(-3,3,,100)
mle <- coef(fm.occ)
bayesB <- out1$summary[1:4,1]
bayesJ <- out2$BUGSoutput$summary[1:4,1]
plot(vegHt, z, xlab="Vegetation height", ylab="", main="Occupancy probability", las = 1, lty = 1)
lines(vegHt, psi, lwd=3, col="black")
lines(predHt, plogis(mle[1] + mle[2]*predHt), lwd = 3, col="blue", lty = 1)
lines(predHt, plogis(bayesB[1] + bayesB[2]*predHt), lwd = 3, col="green", lty = 2)
lines(predHt, plogis(bayesJ[1] + bayesJ[3]*predHt), lwd = 3, col="red", lty = 3)

plot(wind, y, xlab="Wind", ylab="", main="Detection probability", las = 1, lty = 1)
lines(sort(wind), p[order(wind)], lwd=3, col="black")
lines(predW, plogis(mle[3] + mle[4]*predW), lwd = 3, col="blue", lty = 1)
lines(predW, plogis(bayesB[3] + bayesB[4]*predW), lwd = 3, col="green", lty = 2)
lines(predW, plogis(bayesJ[2] + bayesJ[4]*predW), lwd = 3, col="red", lty = 3)
legend(-2.5, 0.3, c("Truth", "unmarked", "WinBUGS", "JAGS"), col=c("black", "blue", "green", "red"), lty=c(1, 1, 2, 3), lwd=3)


### 8.5. Analysis of the distribution of red squirrels in Switzerland using unmarked
# Read in data set
cs <- read.table("crossbill.squirrel.txt", header = TRUE)
str(cs)

# Grab 2006 data for analysis
y <- as.matrix(cs[cs$spec.name == "Squirrel", 7:9])
str(y)
ele <- cs[cs$spec.name == "Squirrel", "ele"]
forest <- cs[cs$spec.name == "Squirrel", "forest"]
date <- as.matrix(cs[cs$spec.name == "Squirrel", 13:15])
dur <- as.matrix(cs[cs$spec.name == "Squirrel", 19:21])

# Standardize covariates and mean-impute date
ele.mean <- mean(ele)
ele.sd <- sd(ele)
forest.mean <- mean(forest)
forest.sd <- sd(forest)
date.mean <- mean(c(date), na.rm = TRUE)
date.sd <- sd(c(date), na.rm = TRUE)
ele <- (ele - ele.mean)/ele.sd
forest <- (forest - forest.mean)/forest.sd
date <- (date - date.mean) / date.sd
date[is.na(date)] <- 0

# Load library, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(ele = ele, forest = forest), obsCovs = list(date = date))
summary(umf)

# Fit a series of models for detection first and do model selection
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~date ~1, data=umf))
summary(fm3 <- occu(~date+I(date^2) ~1, data=umf))
summary(fm4 <- occu(~date+I(date^2)+I(date^3) ~1, data=umf))

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(.)"                    = fm1,
               "p(date)psi(.)"                 = fm2,
               "p(date+date2)psi(.)"           = fm3,
               "p(date+date2+date3)psi(.)"     = fm4)

(ms <- modSel(fms))

# Continue with model fitting for occupancy, guided by AIC as we go
# Effects of elevation
summary(fm5 <- occu(~date+I(date^2)+I(date^3) ~ele, data=umf))
summary(fm6 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2), data=umf))
summary(fm7 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+ I(ele^3), data=umf))
# ele2 best

# Effects of forest and interactions
summary(fm8 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest, data=umf))
summary(fm9 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2), data=umf))
summary(fm10 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2)+ele:forest, data=umf))
summary(fm11 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2)+ele:forest+ele:I(forest^2), data=umf))
summary(fm12 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2)+ele:forest+ele:I(forest^2)+I(ele^2):forest, data=umf))
summary(fm13 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2)+ele:forest+ele:I(forest^2)+I(ele^2):forest+ I(ele^2):I(forest^2), data=umf))
summary(fm14 <- occu(~date+I(date^2)+I(date^3) ~ele+I(ele^2)+forest+I(forest^2)+ele:forest+ele:I(forest^2)+I(ele^2):forest+ I(ele^2):I(forest^2), data=umf))
# fm8 best and most sensible biologically

# Check for additional effects in detection
summary(fm15 <- occu(~date+I(date^2)+I(date^3)+ele ~ele+I(ele^2)+forest, data=umf))
summary(fm16 <- occu(~date+I(date^2)+I(date^3)+ele+I(ele^2) ~ele+I(ele^2)+forest, data=umf))
# ele sufficient in p, check for forest
summary(fm17 <- occu(~date+I(date^2)+I(date^3)+ele+forest ~ele+I(ele^2)+forest, data=umf))
# No effect of forest on p

# Create new covariates for prediction and get predictions
orig.ele <- seq(200, 2500,,100)
orig.forest <- 0:100
orig.date <- seq(15, 110,,100)

ele.pred <- (orig.ele - ele.mean) / ele.sd
forest.pred <- (orig.forest - forest.mean) / forest.sd
date.pred <- (orig.date - date.mean) / date.sd

beta <- coef(fm15)
pred.occ.ele <- plogis(beta[1] + beta[2]*ele.pred + beta[3]*ele.pred^2)
pred.occ.forest <- plogis(beta[1] + beta[4]*forest.pred)
pred.det.date <- plogis(beta[5] + beta[6]*date.pred + beta[7]*date.pred^2 + beta[8]*date.pred^3)
pred.det.ele <- plogis(beta[5] + beta[9]*ele.pred)


par(mfrow = c(2,2), mar = c(5,5,2,3), cex.lab = 1.2)
plot(pred.occ.ele ~ orig.ele, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Predicted occupancy probability", xlab = "Elevation (m)")
plot(pred.occ.forest ~ orig.forest, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Predicted occupancy probability", xlab = "Forest cover (%)")
plot(pred.det.date ~ orig.date, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Predicted detection probability", xlab = "Date (1 = 1 April)")
plot(pred.det.ele ~ orig.ele, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Predicted detection probability", xlab = "Elevation (m)")

# Get the Swiss landscape data again
landscape <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/Swiss_landscape.csv?attredirects=0&d=1")
# landscape<-read.csv("Swiss_landscape.csv")
head(landscape)

# ------------------ Do 2D predictions (maps) ----------------------------#
library(sp)
# Get coordinates and rescale from km to m
coordCH <- matrix(cbind(landscape$x + 0.5, landscape$y + 0.5), ncol = 2)
xcor <- coordCH[,1] * 1000
ycor <- coordCH[,2] * 1000

# Get predictions of occupancy for each 1km2 quadrat of Switzerland
# First scale values of elevation analogously to analysis
pelev <- (landscape$medel-ele.mean)/ele.sd
pforest <- (landscape$forest-forest.mean)/forest.sd

par(mfrow=c(1,1))
betavec<-beta[1:4]               # Coefficients in occupancy
Xg<-cbind(rep(1,length(pelev)),pelev,pelev^2, pforest)
pred<-plogis(Xg%*%(betavec))

# Define a new dataframe with coordinates and outcome to be plotted 
PARAM <- data.frame(x = xcor, y = ycor, z = pred)

# Convert the dataframe first into a SpatialPointsDataFrame and then into a SpatialPixelsDataFrame
coordinates(PARAM)<- ~x+y
gridded(PARAM) <- TRUE

# Plot the map using custom color palette
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
spplot(PARAM, col.regions = mapPalette(100), main = "Red squirrel occupancy probability 2006")

# Get area of squirrel occurrence in 2006
sum(pred)


# Get predictions of detection for each 1km2 quadrat of Switzerland
betavec<-beta[c(5,9)]               # Coefficient of elevation in detection
Xg<-cbind(rep(1,length(pelev)), pelev)
pred<-plogis(Xg%*%(betavec))

# Define a new dataframe with coordinates and outcome to be plotted 
PARAM <- data.frame(x = xcor, y = ycor, z = pred)

# Convert the dataframe first into a SpatialPointsDataFrame and then into a SpatialPixelsDataFrame
coordinates(PARAM)<- ~x+y
gridded(PARAM) <- TRUE

# Plot the map using custom color palette
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
spplot(PARAM, col.regions = mapPalette(100), main = "Red squirrel detection probability 2006")


### 8.6. Two-species occupancy models with asymmetric interaction in BUGS/JAGS

# Grab and bundle data
sq <- as.matrix(cs[cs$spec.name == "Squirrel", 7:9])
cb <- as.matrix(cs[cs$spec.name == "Crossbill", 7:9])
win.data <- list(sq = sq, cb = cb, R = nrow(sq), T = ncol(sq))


# Specify model in BUGS language
sink("model.txt")
cat("
model {

# ----------------- Model for squirrels ------------------
# Priors
psi.sq ~ dunif(0,1)
p.sq ~ dunif(0,1)

# Likelihood
for (i in 1:R) {
   z.sq[i] ~ dbern(psi.sq)
   for (j in 1:T) {
      sq[i,j] ~ dbern(p.eff.sq[i,j])
      p.eff.sq[i,j] <- z.sq[i] * p.sq
   }
}
occ.fs.sq <- sum(z.sq[])

# ----------------- Model for crossbills ------------------
# Priors
psi.cb.with ~ dunif(0,1)
psi.cb.without ~ dunif(0,1)
p.cb.with ~ dunif(0,1)
p.cb.without ~ dunif(0,1)

# Likelihood
for (i in 1:R) {
   z.cb[i] ~ dbern(psi.cb[i])
   psi.cb[i] <- z.sq[i] * psi.cb.with + (1 - z.sq[i]) * psi.cb.without
   p.cb[i] <- z.sq[i] * p.cb.with + (1 - z.sq[i]) * p.cb.without
   for (j in 1:T) {
      cb[i,j] ~ dbern(p.eff.cb[i,j])
      p.eff.cb[i,j] <- z.cb[i] * p.cb[i]
   }
}
occ.fs.cb <- sum(z.cb[])
}
",fill = TRUE)
sink()


# Initial values
zst.sq <- apply(sq, 1, max, na.rm = TRUE)
zst.cb <- apply(cb, 1, max, na.rm = TRUE)
zst.sq[zst.sq == -Inf] <- 0
zst.cb[zst.cb == -Inf] <- 0
inits <- function(){list(z.sq = zst.sq, z.cb = zst.cb)}

# Parameters monitored
params <- c("psi.sq", "p.sq", "psi.cb.with", "psi.cb.without", "p.cb.with", "p.cb.without", "occ.fs.sq", "occ.fs.cb")
# add "z.sq" and "z.cb" for random effects estimates

# MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 3

# Call WinBUGS from R (ART 1.0 min)
out1 <- bugs(win.data, inits, params, "model.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 2)

# Call JAGS from R (ART 0.6 min)
library("R2jags")		# requires rjags
system.time(out2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb) )
traceplot(out2)

# Summarize posteriors
print(out2, dig = 2)
