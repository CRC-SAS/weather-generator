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
setwd("/Users/andrew/Desktop/workshop-BA/3_wgen-conditional-seasonal")
## load the setup data
load("../1_setup/wgen-setup.RData")

## compute regional average seasonal total precip
## to use as predictor
montot = array(dim=c(length(uyr),12,np))
for(i in 1:np){
  tmp = as.data.table(cbind(yr,mo,da,PP[,i]))
    for(j in 1:length(unique(mo))){
      tmp1 = tmp[mo==j]
      montot[,j,i] = as.vector(t(ddply(tmp1, .(yr), summarise, total=sum(V4), na.rm=T)[2]))
    }
}

allseatot = array(dim=c(length(uyr),4,np))
for(j in 1:np){
  for(i in 1:length(uyr)){
      allseatot[i,1,j] = sum(montot[i,1:3,j])
      allseatot[i,2,j] = sum(montot[i,4:6,j])
      allseatot[i,3,j] = sum(montot[i,7:9,j])
      allseatot[i,4,j] = sum(montot[i,10:12,j])
  }
}

# calculate averages but do not consider withheld stations
#avgmontot = apply(montot,1:2,mean,na.rm=T)
avgmontot = apply(montot,1:2,mean,na.rm=T)

seatot = matrix(NA,nr=length(uyr),nc=4)
for(i in 1:length(uyr)){
     seatot[i,1] = sum(avgmontot[i,1:3])
     seatot[i,2] = sum(avgmontot[i,4:6])
     seatot[i,3] = sum(avgmontot[i,7:9])
     seatot[i,4] = sum(avgmontot[i,10:12])
}

## create seasonal total vector
## (repeats for every day within season
sealen = c(sum(c(31,28,31)),sum(c(30,31,30)),sum(c(31,31,30)),sum(c(31,30,31)))
sealen.leap = c(sum(c(31,29,31)),sum(c(30,31,30)),sum(c(31,31,30)),sum(c(31,30,31)))
ST = c()
for(i in 1:nrow(seatot)){
  for(j in 1:ncol(seatot)){
    if(yearID[i]==365) ST = c(ST,rep(seatot[i,j],each=sealen[j])) else ST = c(ST,rep(seatot[i,j],each=sealen.leap[j]))
  }
}

ST1 = ST2 = ST3 = ST4 = ST
ST1[mo!=1 & mo!=2 & mo!=3] = 0
ST2[mo!=4 & mo!=5 & mo!=6] = 0
ST3[mo!=7 & mo!=8 & mo!=9] = 0
ST4[mo!=10 & mo!=11 & mo!=12] = 0

# regional mean maximum temperatures as covariates -- by season
maxmean = array(dim=c(length(uyr),12,np))
for(i in 1:np){
  tmp = as.data.table(cbind(yr,mo,da,MX[,i]))
    for(j in 1:length(unique(mo))){
      tmp1 = tmp[mo==j]
      maxmean[,j,i] = as.vector(t(ddply(tmp1, .(yr), summarise, total=mean(V4), na.rm=T)[2]))
    }
}

allseamax = array(dim=c(length(uyr),4,np))
for(j in 1:np){
  for(i in 1:length(uyr)){
      allseamax[i,1,j] = mean(maxmean[i,1:3,j])
      allseamax[i,2,j] = mean(maxmean[i,4:6,j])
      allseamax[i,3,j] = mean(maxmean[i,7:9,j])
      allseamax[i,4,j] = mean(maxmean[i,10:12,j])
  }
}

# calculate averages but do not consider withheld stations
avgmaxmean = apply(maxmean,1:2,mean,na.rm=T)

seamax = matrix(NA,nr=length(uyr),nc=4)
for(i in 1:length(uyr)){
     seamax[i,1] = mean(avgmaxmean[i,1:3])
     seamax[i,2] = mean(avgmaxmean[i,4:6])
     seamax[i,3] = mean(avgmaxmean[i,7:9])
     seamax[i,4] = mean(avgmaxmean[i,10:12])
}

## create seasonal total vector (repeats for every day within season)
SMX = c()
for(i in 1:nrow(seamax)){
  for(j in 1:ncol(seamax)){
    if(yearID[i]==365) SMX = c(SMX,rep(seamax[i,j],each=sealen[j])) else SMX = c(SMX,rep(seamax[i,j],each=sealen.leap[j]))
  }
}

SMX1 = SMX2 = SMX3 = SMX4 = SMX
SMX1[mo!=1 & mo!=2 & mo!=3] = 0
SMX2[mo!=4 & mo!=5 & mo!=6] = 0
SMX3[mo!=7 & mo!=8 & mo!=9] = 0
SMX4[mo!=10 & mo!=11 & mo!=12] = 0

# regional mean minimum temperatures as covariates -- by season
minmean = array(dim=c(length(uyr),12,np))
for(i in 1:np){
  tmp = as.data.table(cbind(yr,mo,da,MN[,i]))
    for(j in 1:length(unique(mo))){
      tmp1 = tmp[mo==j]
      minmean[,j,i] = as.vector(t(ddply(tmp1, .(yr), summarise, total=mean(V4), na.rm=T)[2]))
    }
}

allseamin = array(dim=c(length(uyr),4,np))
for(j in 1:np){
  for(i in 1:length(uyr)){
      allseamin[i,1,j] = mean(minmean[i,1:3,j])
      allseamin[i,2,j] = mean(minmean[i,4:6,j])
      allseamin[i,3,j] = mean(minmean[i,7:9,j])
      allseamin[i,4,j] = mean(minmean[i,10:12,j])
  }
}

# calculate averages but do not consider withheld stations
avgminmean = apply(minmean,1:2,mean,na.rm=T)

seamin = matrix(NA,nr=length(uyr),nc=4)
for(i in 1:length(uyr)){
     seamin[i,1] = mean(avgminmean[i,1:3])
     seamin[i,2] = mean(avgminmean[i,4:6])
     seamin[i,3] = mean(avgminmean[i,7:9])
     seamin[i,4] = mean(avgminmean[i,10:12])
}

## create seasonal total vector (repeats for every day within season)
SMN = c()
for(i in 1:nrow(seamin)){
  for(j in 1:ncol(seamin)){
    if(yearID[i]==365) SMN = c(SMN,rep(seamin[i,j],each=sealen[j])) else SMN = c(SMN,rep(seamin[i,j],each=sealen.leap[j]))
  }
}

SMN1 = SMN2 = SMN3 = SMN4 = SMN
SMN1[mo!=1 & mo!=2 & mo!=3] = 0
SMN2[mo!=4 & mo!=5 & mo!=6] = 0
SMN3[mo!=7 & mo!=8 & mo!=9] = 0
SMN4[mo!=10 & mo!=11 & mo!=12] = 0


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
	X.OCC <- cbind(POCC[,i],ct,st,ST1,ST2,ST3,ST4)
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
	X.AMT <- cbind(ct[missing.id], st[missing.id],ST1[missing.id],ST2[missing.id],ST3[missing.id],ST4[missing.id])
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
	X.MIN <- cbind(PMN[,i], PMX[,i], ct, st, OCC[,i], Rt, SMN1, SMN2, SMN3, SMN4, SMX1, SMX2, SMX3, SMX4)
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
	X.MAX <- cbind(PMN[,i], PMX[,i], ct, st, OCC[,i], Rt, SMN1, SMN2, SMN3, SMN4, SMX1, SMX2, SMX3, SMX4)
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

# we calculate the covariance of maximum and minimum temperature residuals
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
## in this tutorial we simulate an arbitrary number of trajectories of OND 2015
## 
OND.length <- 31+30+31
yr.sim <- rep(2015,OND.length)
mo.sim <- c(rep(10,31),rep(11,30),rep(12,31))
da.sim <- c(1:31,1:30,1:31)
ct.sim <- ct[(nt-OND.length+1):nt]
st.sim <- st[(nt-OND.length+1):nt]
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
NT=100
RNum <- 1:NT
n.realizations <- length(RNum)
# id for missing data
mv <- -9999

# bootstrapping OND (i.e., ST4, SMN4, SMX4) values to condition output
# OND 2015 FORECAST FOR PRECIPITATION: 70:20:10
# OND 2015 FORECAST FOR TEMPERATURE:   40:35:25
# see corresponding gif files for IRI forecasts

# sample NT different OND precip totals with IRI probability (B, N, A) = (0.1, 0.2, 0.7)
st.samp = sample(1:3, NT, prob=c(0.1, 0.2, 0.70), replace=TRUE)

st.levels = list()
st.levels[[1]] = seatot[,4][seatot[,4] <= quantile(seatot[,4],(1/3))]
st.levels[[2]] = seatot[,4][seatot[,4] <= quantile(seatot[,4],(2/3)) & seatot[,4] > quantile(seatot[,4],(1/3))]
st.levels[[3]] = seatot[,4][seatot[,4] >  quantile(seatot[,4],(2/3))]

# we are only simulating ST4 (OND), so ST1 = ST2 = ST3 = rep(0, season.length) ... where season.length differs between ST1, ST2, and ST3
ST1.sim = ST2.sim = ST3.sim = ST4.sim = matrix(0,nrow=NT,ncol=length(ct.sim))
for(k in 1:NT) ST4.sim[k,] = rep(sample(st.levels[[st.samp[k]]],1),length(ct.sim))

# sample NT different OND max and min temp totals with IRI probability (B, N, A) = (0.25, 0.35, 0.40)
tmp.samp = sample(1:3, NT, prob=c(0.25, 0.35, 0.40), replace=TRUE)
mx.levels = list()
mx.levels[[1]] = seamax[,4][seamax[,4] <= quantile(seamax[,4],(1/3))]
mx.levels[[2]] = seamax[,4][seamax[,4] <= quantile(seamax[,4],(2/3)) & seamax[,4] > quantile(seamax[,4],(1/3))]
mx.levels[[3]] = seamax[,4][seamax[,4] >  quantile(seamax[,4],(2/3))]

# we are only simulating SMX4 (OND), so SMX1 = SMX2 = SMX3 = rep(0, season.length) ... where season.length differs between SMX1, SMX2, and SMX3
SMX1.sim = SMX2.sim = SMX3.sim = SMX4.sim = matrix(0,nrow=NT,ncol=length(ct.sim))
for(k in 1:NT) SMX4.sim[k,] = rep(sample(mx.levels[[tmp.samp[k]]],1),length(ct.sim))

mn.levels = list()
mn.levels[[1]] = seamin[,4][seamin[,4] <= quantile(seamin[,4],(1/3))]
mn.levels[[2]] = seamin[,4][seamin[,4] <= quantile(seamin[,4],(2/3)) & seamin[,4] > quantile(seamin[,4],(1/3))]
mn.levels[[3]] = seamin[,4][seamin[,4] >  quantile(seamin[,4],(2/3))]

# we are only simulating SMN4 (OND), so SMN1 = SMN2 = SMN3 = rep(0, season.length) ... where season.length differs between SMN1, SMN2, and SMN3
SMN1.sim = SMN2.sim = SMN3.sim = SMN4.sim = matrix(0,nrow=NT,ncol=length(ct.sim))
for(k in 1:NT) SMN4.sim[k,] = rep(sample(mn.levels[[tmp.samp[k]]],1),length(ct.sim))


# define matrices for simulated series
SIMocc <- SIMamt <- SIMmax <- SIMmin <- array(dim=c(NT,nt.sim,np)) # simulations in array of (# days) x (# stations) x (# trajs)

## GAMMA MODEL -- OBTAIN SHAPE AND SCALE
SC <- matrix(NA,nrow=nt.sim,ncol=np)
SH <- numeric(np)

# can't have max temp less than min temp
#min.diff <- min(MX-MN,na.rm=T)
min.diff <- quantile(MX-MN,probs=1/10000,na.rm=T)
SIMocc.old.2 <- (mvrnorm(n=1,mu=rep(0,np),Sigma=PRCPcor[[mo[nt.sim]]]) > 0) + 0

# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
for(i in 1:NT){
	# this loop is redefining the Gamma parameters for each sample
	for(k in 1:np){
  		SH[k] <- gamma.shape(GAMMA[[k]])$alpha
  		SC[,k] <- exp(apply(coefamt[,k]*rbind(1,ct.sim,st.sim,ST1.sim[i,],ST2.sim[i,],ST3.sim[i,],ST4.sim[i,]),FUN=sum,MAR=2,na.rm=T))/SH[k]
	}
	SIMocc.old <- SIMocc.old.2
	SIMmax.old <- apply(MX[mo==12 & da==31,],2,mean,na.rm=T)
	SIMmin.old <- apply(MN[mo==12 & da==31,],2,mean,na.rm=T)
	for(d in 1:nt.sim){ 
  		# amounts are transformed from mvrnorm with PRCPcor as covariance matrix 
  		SIMamt[i,d,] <- qgamma(pnorm(mvrnorm(n=1,mu=rep(0,np),Sigma=PRCPcor[[mo.sim[d]]])),shape=SH,scale=SC[d,])
  		# occurrence is mean function + mvrnorm with PRCPcor as covariance matrix (coerced to 0 or 1)
  		# covariates for mean function are POCC, cos[d], sin[d]
  		SIMocc[i,d,] <- (apply(rbind(1,SIMocc.old,ct.sim[d],st.sim[d],ST1.sim[i,d],ST2.sim[i,d],ST3.sim[i,d],ST4.sim[i,d])*coefocc,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=PRCPcor[[mo.sim[d]]]) > 0) + 0
  		# max is mean function + mvrnorm with TEMPcov as covariance matrix
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmax[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmax,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMAXcov[[mo.sim[d]]])
  		# min is mean function + mvrnorm with TEMPcov as covariance matrix (different mvrnorm than max)
  		# covariates for mean function are PMN, PMX, ct, st, OCC, Rt
  		SIMmin[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmin,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMINcov[[mo.sim[d]]])
  		# can't have max temp less than min temp
  		while(min(SIMmax[i,d,] - SIMmin[i,d,]) < min.diff){
  			SIMmin[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmin,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMINcov[[mo.sim[d]]])
  			SIMmax[i,d,] <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc[i,d,],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmax,2,sum) + mvrnorm(1,mu=rep(0,np),Sigma=TMAXcov[[mo.sim[d]]])
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

STs <- ncvar_def(name="ST", units="mm season-1", dim=list(dimRNum), missval=mv, longname="Seasonal total precipitation covariates", prec="float", verbose=TRUE)
SMXs <- ncvar_def(name="SMX", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average maximum temperature covariates", prec="float", verbose=TRUE)
SMNs <- ncvar_def(name="SMN", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average minimum temperature covariates", prec="float", verbose=TRUE)

# Create NetCDF file in current working directory
dd <- nc_create(filename="wgen-conditional-seasonal.nc",
  vars = list(tmax, tmin, prcp, srad, et0, dec_lon, dec_lat, elevation, STs, SMXs, SMNs),
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

ncvar_put(dd, varid = "ST",
  vals = apply(ST4.sim, 1, unique),
  start = 1, count = NT,
  verbose = TRUE)

ncvar_put(dd, varid = "SMX", 
  vals = apply(SMX4.sim, 1, unique),
  start = 1, count = NT,
  verbose = TRUE)

ncvar_put(dd, varid = "SMN",
  vals = apply(SMN4.sim, 1, unique),
  start = 1, count = NT, 
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

ncatt_put( dd, varid = 0,
  attname = "climate driver covariates",
  attval = "Regionally-averaged seasonal total precipitation, average maximum temperature, and average minimum temperature",
  verbose = TRUE)

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
