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
setwd("/Users/andrew/Desktop/workshop-BA/2_wgen-unconditional")
## load the setup data
load("../1_setup/wgen-setup.RData")

## precipitation occurrence is modeled using probit regression
## fit glm at each location separately, save the coefficients
# save the model fit in a list, one list element for each location
PROBIT <- list()
# save the model residuals in a matrix of same dimensions as OCC
PROBITres <- matrix(NA,nrow=nrow(OCC),ncol=ncol(OCC))
# begin loop over locations
for(i in 1:np){
	# define response
	Y.OCC <- OCC[,i]
	# define design matrix (covariates)
	X.OCC <- cbind(POCC[,i],ct,st)
	# save model in list element "i"
	PROBIT[[i]] <- glm(Y.OCC ~ X.OCC, family=binomial(probit))
	# save model residuals in column "i" of PROBITres matrix
	# but there are missing values, so must identify those 
	missing.id <- (!is.na(Y.OCC) & !is.na(apply(X.OCC,1,sum)))
	PROBITres[missing.id,i] <- PROBIT[[i]]$residuals
}
coefocc <- matrix(NA,nrow=length(PROBIT[[1]]$coefficients),ncol=np)
for(i in 1:np){
	coefocc[,i] <- PROBIT[[i]]$coefficients 
}

## precipitation amount is modeled using Gamma model with spatially-varying shape,
## and spatio-temporally varying scale 
# save the model fit in a list, one list element for each location
GAMMA <- list()
# no model residuals, assumed to be negligible
# begin loop over locations
for(i in 1:np){
	# identify which values are NA, remove them (Gamma glm hates NA)
	missing.id <- (!is.na(PPI[,i]))
	# define response
	Y.AMT <- PPI[missing.id,i]
	# define design matrix (covariates)
	X.AMT <- cbind(ct[missing.id], st[missing.id])
	# save model in list element "i"
	GAMMA[[i]] <- glm(Y.AMT ~ X.AMT, family=Gamma(link=log))
}
coefamt <- matrix(NA,nrow=length(GAMMA[[1]]$coefficients),ncol=np)
for(i in 1:np){
	coefamt[,i] <- GAMMA[[i]]$coefficients 
}

## minimum temperature is modeled using linear regression
# save the model fit in a list, one list element for each location
TMIN <- list()
# save the model residuals in a matrix of same dimensions as MN
TMINres <- matrix(NA,nrow=nrow(MN),ncol=ncol(MN))
# begin loop over locations
for(i in 1:np){
	# define response
	Y.MIN <- MN[,i]
	# define design matrix (covariates)
	X.MIN <- cbind(PMN[,i], PMX[,i], ct, st, OCC[,i], Rt)
	# save model in list element "i"
	TMIN[[i]] <- lm(Y.MIN ~ X.MIN)
	# save model residuals in column "i" of TMINres matrix
	# but there are missing values, so must identify those
	missing.id <- (!is.na(Y.MIN) & !is.na(apply(X.MIN,1,sum)))
	TMINres[missing.id,i] <- TMIN[[i]]$residuals 
}
coefmin <- matrix(NA,nrow=length(TMIN[[1]]$coefficients),ncol=np)
for(i in 1:np){
	coefmin[,i] <- TMIN[[i]]$coefficients 
}

## maximum temperature is odeled using linear regression
# save the model fit in a list, one list element for each location
TMAX <- list()
# save the model residuals in a matrix of same dimensions as MX
TMAXres <- matrix(NA,nrow=nrow(MX),ncol=ncol(MX))
# begin loop over locations
for(i in 1:np){
	# define response
	Y.MAX <- MX[,i]
	# define design matrix (covariates)
	X.MAX <- cbind(PMN[,i], PMX[,i], ct, st, OCC[,i], Rt)
	# save model in list element "i"
	TMAX[[i]] <- lm(Y.MAX ~ X.MAX)
	# save model residuals of column "i" of TMAXres matrix
	# but there are missing values, so must identify those
	missing.id <- (!is.na(Y.MAX) & !is.na(apply(X.MAX,1,sum)))
	TMAXres[missing.id,i] <- TMAX[[i]]$residuals 
}
coefmax <- matrix(NA,nrow=length(TMAX[[i]]$coefficients),ncol=np)
for(i in 1:np){
	coefmax[,i] <- TMAX[[i]]$coefficients 
}

# estimate the spatial covariance of residuals for maximum and minimum temperatures
TMAXcov <- TMINcov <- list()
# begin loop over months, to account for seasonality trends in model residuals
for(i in 1:12){
	TMAXcov[[i]] <- cov(TMAXres[mo==i,],use="complete")
  TMINcov[[i]] <- cov(TMINres[mo==i,],use="complete")
}
# we now much calculate the spatial correlation matrix for precipitation
# occurrence residuals...  CORRELATION matrices are used instead of COVARIANCE
# matrices because probit regression has variance unity by definition.
PRCPcor <- list()
# begin loop over months, to account for seasonality trends in model residuals
for(i in 1:12){
	PRCPcor[[i]] <- cor(PROBITres[mo==i,],use="pairwise.complete")
}


##########################################
##
## Simulations
##
##########################################
# simulation metadata
x_coords = lon.lat[,1]
y_coords = lon.lat[,2]
n.x_coords <- length(x_coords)
n.y_coords <- length(y_coords)
# elevation at each station 
elev <- stns[,"elev"]
## in this tutorial we simulate an arbitrary number of trajectories of arbitrary length
## let us assume the two years we are simulating are NOT leap years...2017-2018
## 
## now... the first two years in the record (from setup file) are 1961-1962,
## both of which are not leap years.  therefore, we can use the first two 
## years to define our simulation months and days, and subsequently sine and cosine.
##  we define years as 2017-2018
yr.sim <- rep(2017:2018,each=365)
mo.sim <- mo[1:length(yr.sim)]
da.sim <- da[1:length(yr.sim)]
ct.sim <- ct[1:length(yr.sim)]
st.sim <- st[1:length(yr.sim)]
## if you wish to neglect temperature trends, and be centered on the mean of 1961-2013 
## temperatures, set Rt.sim = 0 for all days
Rt.sim <- rep(0,length(yr.sim))
#Rt.sim <- rep(1,length(yr.sim)) # if you want simulated temperatures to be high
#Rt.sim <- rep(-1,length(yr.sim)) # if you want simulated temperatures to be low

# number of days to simulate
nt.sim <- length(Rt.sim)

## years to simulate (just for convenience, they are in fact arbitrary years)
uyr.sim <- unique(yr.sim)
## number of years to simulate (also arbitrary)
nyr.sim <- length(uyr.sim)

# simulating 1 Jan 2017 -- 31 Dec 2018
sim.start.date <- as.Date(paste(yr.sim[1],"-",mo.sim[1],"-",da.sim[1],sep=""))   # Beginning of simulated series
sim.end.date <- as.Date(paste(yr.sim[length(yr.sim)],"-",mo.sim[length(mo.sim)],"-",da.sim[length(da.sim)],sep=""))     # End of simulated series
sim.dates <- seq(from = sim.start.date, to = sim.end.date, by = "days")
# julian days (for .nc file)
Times <- julian(sim.dates, origin = as.Date("1961-01-01"))
n.times <- length(Times)    # Number of simulated days
# number of realizations 
NT=10#0
RNum <- 1:NT
n.realizations <- length(RNum)
# id for missing data
mv <- -9999

# define matrices for simulated series
SIMocc <- SIMamt <- SIMmax <- SIMmin <- array(dim=c(NT,nt.sim,np)) # simulations in array of (# days) x (# stations) x (# trajs)

## GAMMA MODEL -- OBTAIN SHAPE AND SCALE
MU <- SC <- matrix(NA,nrow=nt.sim,ncol=np)
SH <- numeric(np)

for(i in 1:np){
  SH[i] <- gamma.shape(GAMMA[[i]])$alpha
  SC[,i] <- exp(apply(GAMMA[[i]]$coef*rbind(1,ct.sim,st.sim),FUN=sum,MAR=2,na.rm=T))/SH[i]
}

# can't have max temp less than min temp
#min.diff <- min(MX-MN,na.rm=T)
min.diff <- quantile(MX-MN,probs=1/10000,na.rm=T)
SIMocc.old.2 <- (mvrnorm(n=1,mu=rep(0,np),Sigma=PRCPcor[[mo[nt.sim]]]) > 0) + 0
# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
for(i in 1:NT){
	SIMocc.old <- SIMocc.old.2
	SIMmax.old <- apply(MX[mo==12 & da==31,],2,mean,na.rm=T)
	SIMmin.old <- apply(MN[mo==12 & da==31,],2,mean,na.rm=T)
	for(d in 1:nt.sim){ 
  		# amounts are transformed from mvrnorm with PRCPcor as covariance matrix 
  		SIMamt[i,d,] <- qgamma(pnorm(mvrnorm(n=1,mu=rep(0,np),Sigma=PRCPcor[[mo.sim[d]]])),shape=SH,scale=SC[d,])
  		# occurrence is mean function + mvrnorm with PRCPcor as covariance matrix (coerced to 0 or 1)
  		# covariates for mean function are POCC, cos[d], sin[d]
  		SIMocc[i,d,] <- (apply(rbind(1,SIMocc.old,ct.sim[d],st.sim[d])*coefocc,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=PRCPcor[[mo.sim[d]]]) > 0) + 0
  		# max is mean function + mvrnorm with TEMPcov as covariance matrix
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmax[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d])*coefmax,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMAXcov[[mo.sim[d]]])
  		# min is mean function + mvrnorm with TEMPcov as covariance matrix (different mvrnorm than max)
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmin[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d])*coefmin,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMINcov[[mo.sim[d]]])
  		# can't have max temp less than min temp
  		while(min(SIMmax[i,d,] - SIMmin[i,d,]) < min.diff){
  			SIMmin[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d])*coefmin,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMINcov[[mo.sim[d]]])
  			SIMmax[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d])*coefmax,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMAXcov[[mo.sim[d]]])
  		}
  		SIMocc.old <- SIMocc[i,d,]
  		SIMmax.old <- SIMmax[i,d,]
  		SIMmin.old <- SIMmin[i,d,]
	}
setTxtProgressBar(pb, i)
}

SIMamt[SIMamt < 0.1] = 0
SIMocc[SIMamt ==  0] = 0
SIMamt[SIMocc ==  0] = 0
 
# estimate srad, et0
# set up for calculating reference evapotranspiration (ET)
  coef.a <- 0.0019  
  coef.b <- 21.0584
  coef.c <- 0.0874
  coef.d <- 0.6278
# set up for calculating solar radiation (SRAD)
  bc.coef=c(A=0.69, B=0.02, C=2.12)
# set up for both ET and SRAD
  stns.rad = radians(lon.lat[,2])
  day.of.year = yday(sim.dates)
  q0 = matrix(NA,nrow=length(day.of.year),ncol=length(stns.rad))
    for(kk in 1:length(stns.rad)) q0[,kk] <- extrat(i=day.of.year, lat=stns.rad[kk])$"ExtraTerrestrialSolarRadiationDaily"

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
extraT <- array(suppressWarnings(extrat(i=matrix(dayOfYear(matrix(rep(sim.dates,np),nrow=nt.sim,ncol=np,byrow=F)),nrow=nt.sim,ncol=np,byrow=F), lat=stns.rad)$ExtraTerrestrialSolarRadiationDaily),dim=dim(dtr))
SIMsrad <- extraT * bc.coef[['A']] * (1 - exp(-bc.coef[['B']] * (dtr^bc.coef[['C']])))

# define dimensions for .nc file 
dimRNum <- ncdim_def(name="rnum", units="number", vals=RNum, unlim=TRUE, create_dimvar=TRUE, longname="Realization number")
dimStation <- ncdim_def(name="stn_id", units="number", vals=1:np, unlim=FALSE, create_dimvar=TRUE, longname="Station number corresponding to metadata")
dimTime <- ncdim_def(name="time", units="days since 1961-01-01", vals=Times, unlim = FALSE, create_dimvar=TRUE, calendar="standard", longname="Time in days since 1961-01-01")

tmax <- ncvar_def(name="tmax", units="degrees_Celsius", dim=list(dimRNum, dimTime, dimStation), missval=mv, longname="Daily maximum near-surface air temperature", prec="float", verbose=TRUE)
tmin <- ncvar_def(name="tmin", units="degrees_Celsius", dim=list(dimRNum, dimTime, dimStation), missval=mv, longname="Daily minimum near-surface air temperature", prec="float", verbose=TRUE)
prcp <- ncvar_def(name="prcp", units="mm day-1", dim=list(dimRNum, dimTime, dimStation), missval=mv, longname="Daily total rainfall", prec="float", verbose=TRUE)
srad <- ncvar_def(name="srad", units="Mjoules m-2 day-1", dim=list(dimRNum, dimTime, dimStation), missval=mv, longname="Daily solar radiation", prec="float", verbose=TRUE)
et0  <- ncvar_def(name="et0",  units="mm day-1", dim=list(dimRNum, dimTime, dimStation), missval=mv, longname="Reference evapotranspiration", prec="float", verbose=TRUE)
dec_lon <- ncvar_def(name="dec_lon", units="degrees_east", missval=mv, dim=list(dimStation), longname="Longitude in decimal degrees", prec="float", verbose=TRUE)
dec_lat <- ncvar_def(name="dec_lat", units="degrees_north", missval=mv, dim=list(dimStation), longname="Latitude in decimal degrees", prec="float", verbose=TRUE)
elevation <- ncvar_def(name="elevation", units="meters from sea level", missval=mv, dim=list(dimStation), longname="Elevation in meters from sea level", prec="float", verbose=TRUE)

# Create NetCDF file in current working directory
dd <- nc_create(filename="wgen-unconditional.nc",
  vars = list(tmax, tmin, prcp, srad, et0, dec_lon, dec_lat, elevation),
  force_v4 = TRUE, verbose = TRUE)

ncvar_put(dd, varid = "tmax", 
  vals = SIMmax,
  start = c(1, 1, 1), count = c(NT, nt.sim, np),
  verbose = TRUE)

ncvar_put(dd, varid = "tmin", 
  vals = SIMmin,
  start = c(1, 1, 1), count = c(NT, nt.sim, np),
  verbose = TRUE)

ncvar_put(dd, varid = "prcp", 
  vals = SIMamt,
  start = c(1, 1, 1), count = c(NT, nt.sim, np),
  verbose = TRUE)
  
ncvar_put(dd, varid = "et0", 
  vals = SIMet0,
  start = c(1, 1, 1), count = c(NT, nt.sim, np),
  verbose = TRUE)

ncvar_put(dd, varid = "srad",
  vals = SIMsrad,
  start = c(1, 1, 1), count = c(NT, nt.sim, np),
  verbose = TRUE)

ncvar_put(dd, varid = "dec_lon",
  vals = stns[,"lon_dec"],
  start = 1, count = np,
  verbose = TRUE)

ncvar_put(dd, varid = "dec_lat",
  vals = stns[,"lat_dec"],
  start = 1, count = np,
  verbose = TRUE)

ncvar_put(dd, varid = "elevation",
  vals = stns[,"elev"],
  start = 1, count = np,
  verbose = TRUE)

# --- Write global attributes

ncatt_put( dd, varid = 0,
  attname = "title",
  attval = "Synthetic weather data for the Salado River Basin A (Argentina)",
  verbose = TRUE)

ncatt_put( dd, varid = 0,
  attname = "software",
  attval = "Stochastic weather generator version 3.0",
  verbose = TRUE)

#ncatt_put( dd, varid = 0,
#  attname = "climate driver covariates",
#  attval = "Regionally-averaged seasonal total precipitation, average maximum temperature, and average minimum temperature",
#  verbose = TRUE)

ncatt_put(dd, varid = 0,
  attname = "Start and end dates",
  attval = paste(sim.start.date, sim.end.date),
  verbose = TRUE)
  
ncatt_put(dd, varid = "et0",
  attname = "calculation",
  attval = "Hargreaves-Samani modified by Droogers and Allen",
  verbose = TRUE)

ncatt_put(dd, varid = "et0",
  attname = "coefficients",
  attval = "Coefficients for daily data calibrated for Junín",
  verbose = TRUE)  

ncatt_put(dd, varid = "srad",
          attname = "calculation",
          attval = "Bristow-Campbell, not considering tomorrow's min temp because obs data derived from hourly aggregate Met Service data",
          verbose = TRUE)

ncatt_put(dd, varid = "srad",
          attname = "coefficients",
          attval = "Coefficients estimated for Buenos Aires and Pilar",
          verbose = TRUE)

nc_close(dd)
