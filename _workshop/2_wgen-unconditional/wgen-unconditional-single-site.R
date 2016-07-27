# --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://lib.stat.cmu.edu/R/CRAN/"), error=traceback)

list.of.packages <- c("ncdf4","MASS","fields","geoR","data.table","plyr","lubridate","sirad")  

for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack)
    require(pack, character.only = TRUE)
  }
}

rm(list.of.packages, pack)

## set the working directory -- this is wherever you put the folder from the zip file
setwd("/home/federico/Desarrollo/RProjects/Wth Simulation/")
## load the setup data
load("../1_setup/wgen-setup.RData")


# let us isolate the station we wish to model
stn.id <- 8 # junin, as can be seen in the object named 'stns'
# precip intensity at junin
PPI <- PPI[,stn.id]
# precip (with zeros) at junin
PP <- PP[,stn.id]
# precip occurrence at junin
OCC <- OCC[,stn.id]
# prev. day's precip occurrence at junin
POCC <- POCC[,stn.id]
# min temp at junin
MN <- MN[,stn.id]
# prev. day's min temp at junin
PMN <- PMN[,stn.id]
# max temp at junin
MX <- MX[,stn.id]
# prev. day's max temp at junin
PMX <- PMX[,stn.id]

# longitude and latitude of location
lon.lat <- lon.lat[stn.id,]


## precipitation occurrence is modeled using probit regression

# define design matrix (covariates)
X.OCC <- cbind(POCC,ct,st)

# save model in object named "PROBIT"
PROBIT <- glm(OCC ~ X.OCC, family=binomial(probit))
	
# save the model residuals in a vector of length equal to observed days
PROBITres <- numeric(length(OCC))

# but there are missing values, so must identify those 
missing.id <- (!is.na(OCC) & !is.na(apply(X.OCC,1,sum)))
PROBITres[missing.id] <- PROBIT$residuals

# save the coefficients
coefocc <- PROBIT$coefficients 


## precipitation amount is modeled using Gamma model with spatially-varying shape,
## and spatio-temporally varying scale 

# no model residuals, assumed to be negligible

# identify which values are NA, remove them (Gamma glm hates NA)
missing.id <- (!is.na(PPI))

# define response
Y.AMT <- PPI[missing.id]

# define design matrix (covariates)
X.AMT <- cbind(ct[missing.id], st[missing.id])
	
# save model in object named "GAMMA"
GAMMA <- glm(Y.AMT ~ X.AMT, family=Gamma(link=log))

# save the coefficients
coefamt <- GAMMA$coefficients 

# save the shape and scale
SH <- gamma.shape(GAMMA)$alpha
SC <- exp(apply(coefamt*rbind(1,ct,st),FUN=sum,MAR=2,na.rm=T))/SH


## minimum temperature is modeled using linear regression

# define design matrix (covariates)
X.MIN <- cbind(PMN, PMX, ct, st, OCC, Rt)

# save model in objecdt named "TMIN"
TMIN <- lm(MN ~ X.MIN)

# save the model residuals in a vector of length equal to observed days
TMINres <- numeric(length(MN))

# but there are missing values, so must identify those
missing.id <- (!is.na(MN) & !is.na(apply(X.MIN,1,sum)))
TMINres[missing.id] <- TMIN$residuals 

# save the coefficients
coefmin <- TMIN$coefficients 


## maximum temperature is modeled using linear regression

# define design matrix (covariates)
X.MAX <- cbind(PMN, PMX, ct, st, OCC, Rt)

# save model in object named "TMAX"
TMAX <- lm(MX ~ X.MAX)

# save the model residuals in a vector of length equal to observed days
TMAXres <- numeric(length(MX))

# but there are missing values, so must identify those
missing.id <- (!is.na(MX) & !is.na(apply(X.MAX,1,sum)))
TMAXres[missing.id] <- TMAX$residuals 

# save the coefficients
coefmax <- TMAX$coefficients 


# we want standard deviations of the residuals for temperatures for each month
TMAX.sd <- TMIN.sd <- numeric(12)
for(i in 1:12){
  TMAX.sd[i] <- sd(TMAXres[mo==i],na.rm=T)
  TMIN.sd[i] <- sd(TMINres[mo==i],na.rm=T)
}


##########################################
##
## Simulations
##
##########################################
# number of realizations 
NT=100

# let's simulate the first five years of the observational period (1 Jan 1961 - 31 Dec 1965)
# these are already defined from the setup script
yr.sim <- yr[1:which(yr==1965 & mo==12 & da==31)]
mo.sim <- mo[1:which(yr==1965 & mo==12 & da==31)]
da.sim <- da[1:which(yr==1965 & mo==12 & da==31)]
ct.sim <- ct[1:which(yr==1965 & mo==12 & da==31)]
st.sim <- st[1:which(yr==1965 & mo==12 & da==31)]
Rt.sim <- Rt[1:which(yr==1965 & mo==12 & da==31)]

# number of days to simulate
nt.sim <- length(Rt.sim)

## years to simulate (just for convenience, they are in fact arbitrary years)
uyr.sim <- unique(yr.sim)
## number of years we simulate
nyr.sim <- length(uyr.sim)

# simulating 1 Jan 1961 -- 31 Dec 1965
sim.start.date <- as.Date(paste(yr.sim[1],"-",mo.sim[1],"-",da.sim[1],sep=""))   # Beginning of simulated series
sim.end.date <- as.Date(paste(yr.sim[length(yr.sim)],"-",mo.sim[length(mo.sim)],"-",da.sim[length(da.sim)],sep=""))     # End of simulated series
sim.dates <- seq(from = sim.start.date, to = sim.end.date, by = "days")
# julian days (for .nc file)
Times <- julian(sim.dates, origin = as.Date("1961-01-01"))
n.times <- length(Times)    # Number of simulated days

# id for missing data
mv <- -9999

# define matrices for simulated series
SIMocc <- SIMamt <- SIMmax <- SIMmin <- array(dim=c(nt.sim,NT)) # simulations in array of (# days x # trajectories)


# can't have max temp less than min temp
min.diff <- quantile(MX-MN,probs=1/10000,na.rm=T)

# set the progress bar to know how long this will take
pb <- txtProgressBar(min = 0, max = NT, style = 3)

# being loop over realizations
for(i in 1:NT){

  # initialize with a sample from binomial with probability equal to 
  # the observed probability of occurrence for 31 Dec (all years)
  SIMocc.old <- rbinom(n=1, size=1, prob=mean(OCC[mo==12 & da==31],na.rm=T))

  # initialize temperatures with the climatological mean of 31 Dec
	SIMmax.old <- mean(MX[mo==12 & da==31],na.rm=T)
	SIMmin.old <- mean(MN[mo==12 & da==31],na.rm=T)

	for(d in 1:nt.sim){ 

  		# amounts are transformed from occurrence noise
      # probit has variance 1 by definition, so we use a copula
      # to transform a mean zero, variance unity random sample
  		SIMamt[d,i] <- qgamma(pnorm(rnorm(n=1,mean=0,sd=1)),shape=SH,scale=SC[d])

  		# occurrence is mean function + rnorm(mean=0, sd=1)
      # covariates for mean function are POCC, ct, st
  		SIMocc[d,i] <- ((sum(coefocc*c(1,SIMocc.old,ct.sim[d],st.sim[d])) + rnorm(n=1,mean=0,sd=1)) > 0) + 0

  		# max is mean function + rnorm TMAX.sd according to month
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmax[d,i] <- sum(coefmax*c(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[d,i],Rt.sim[d])) + rnorm(n=1,mean=0,sd=TMAX.sd[mo.sim[d]])

  		# min is mean function + rnorm TMIN.sd according to month
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmin[d,i] <- sum(coefmin*c(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[d,i],Rt.sim[d])) + rnorm(n=1,mean=0,sd=TMIN.sd[mo.sim[d]])

  		# can't have max temp less than min temp
  		while(SIMmax[d,i] - SIMmin[d,i] < min.diff){
  			SIMmin[d,i] <- sum(coefmin*c(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[d,i],Rt.sim[d])) + rnorm(n=1,mean=0,sd=TMIN.sd[mo.sim[d]])
  			SIMmax[d,i] <- sum(coefmax*c(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[d,i],Rt.sim[d])) + rnorm(n=1,mean=0,sd=TMAX.sd[mo.sim[d]])
  		}

      # set previous day covariates as current day's simulated before moving on to next day
  		SIMocc.old <- SIMocc[d,i]
  		SIMmax.old <- SIMmax[d,i]
  		SIMmin.old <- SIMmin[d,i]
	}
setTxtProgressBar(pb, i)
}

SIMamt[SIMamt < 0.1] = 0
SIMocc[SIMamt ==  0] = 0
SIMamt[SIMocc ==  0] = 0
 
# the following coefficients are estimated for the Salado A basin
# you cannot use these for any location!

# estimate srad, et0

# set up for calculating reference evapotranspiration (ET)
  coef.a <- 0.0019  
  coef.b <- 21.0584
  coef.c <- 0.0874
  coef.d <- 0.6278

# set up for calculating solar radiation (SRAD)
  bc.coef=c(A=0.69, B=0.02, C=2.12)

# set up for both ET and SRAD
# convert decimal degrees to radians 
  stns.rad = radians(lon.lat[2])

# define julian days
  day.of.year = yday(sim.dates)

# define vector to hold extra terrestrial solar radiation (length of historic record)
  q0 = numeric(nt.sim)

# begin loop over all days
  for(i in 1:length(q0)){ 
    print(paste0(i,"/",length(q0)))
    q0[i] = as.vector(extrat(i=day.of.year[i], lat=stns.rad)$"ExtraTerrestrialSolarRadiationDaily")[[1]]
  }

# estimate reference ET
# difference between max and min temperatures
dtr <- (SIMmax - SIMmin)
# dtr should NEVER be less than zero (the while statement during simulation ensures this)
dtr[dtr<0] <- NA
# average temperature
tavg = (SIMmax + SIMmin)/2
# equation for reference ET 
SIMet0 <- coef.a * 0.408 * array(q0,dim=dim(dtr)) * (tavg + coef.b) * (dtr - (coef.c * SIMamt)) ** coef.d
# R cant handle imaginary numbers (i.e. exponentiating negative bases)
# so times when temperature gradient is less than modified precipitation amount results in NaN -- set these to 0
SIMet0[is.na(SIMet0)]=0 

# Estimating solar radiation.
# Computate Bristow-Campbell's delta temperature.
# do not consider tomorrow's minimum temperature 
# because this is hourly aggregate Met Service data
# dtr is dtemp
SIMsrad <- array(q0,dim=dim(dtr)) * bc.coef[['A']] * (1 - exp(-bc.coef[['B']] * (dtr^bc.coef[['C']])))

save.image("wgen-unconditional-single-site.RData")
