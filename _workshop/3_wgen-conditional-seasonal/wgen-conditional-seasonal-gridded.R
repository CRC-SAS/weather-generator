# --- Load necessary R packages ----
rm(list=objects()); gc();
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
# setwd("/Users/andrew/Desktop/workshop-BA/3_wgen-conditional-seasonal")
## load the setup data
load("1_setup/wgen-setup.RData")

predlocs <- read.table("1_setup/data/grids.txt",header=T)
predloc.GK <- predlocs[,1:2]/1000 # convert from meters to KILOMETERS
predloc <- predlocs[,3:4]
# distance matrix of station locations
# converts longitude and latitude to KILOMETERS
dist.mat <- rdist.earth(lon.lat,miles=F) 
# diagonal element is not explicitly equal to zero, so define as such
diag(dist.mat) <- 0

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

# we assume residuals from minimum and maximum temperatures come from the
# same distribution, therefore we combine them to estimate the spatial
# covariance matrices...
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

## least squares function for kriging parameters
LS <- function(p){
M <- p[2]*exp((-dist.mat)/p[3])
diag(M) <- p[1]+p[2]
return(sum(vario-M)^2)
}

## estimate kriging parameters
PRCPvario <- TMAXvario <- TMINvario <- list()
params <- params.max <- params.min <- matrix(NA,nrow=12,ncol=3)

for(kk in 1:12){ 
    PRCPvario[[kk]] = var(PROBITres[mo==kk,],use="pairwise.complete")*(1 - PRCPcor[[kk]])
  
  vario <- PRCPvario[[kk]]
  params[kk,] <- optim(par=c(0.01,1,max(dist.mat)),fn=LS)$par

  TMAXvario[[kk]] = cov(TMAXres[mo==kk,],use="pairwise.complete")
  
  vario <- TMAXvario[[kk]]
  params.max[kk,] <- optim(par=c(0.01,mean(var(MX[mo==kk,],na.rm=T)),max(dist.mat)),fn=LS)$par 

  TMINvario[[kk]] = cov(TMINres[mo==kk,],use="pairwise.complete")
  
  vario <- TMINvario[[kk]]
  params.min[kk,] <- optim(par=c(0.01,mean(var(MN[mo==kk,],na.rm=T)),max(dist.mat)),fn=LS)$par
  
  }
params[params<0]=0
params[,1] <- rep(0,12)
params[,2] <- rep(1,12)

params.max[,1] <- params.min[,1] <- rep(0,12)

# estimate model coefficients on grid using ordinary kriging (OK)
# occurrence
coefocc.sim <- matrix(NA,nrow=nrow(coefocc),ncol=nrow(predloc))
for(kk in 1:nrow(coefocc.sim)){
  coefocc.sim[kk,] = suppressWarnings(predict(Krig(lon.lat,coefocc[kk,]),predloc))
}
# minimum temperature
coefmin.sim <- matrix(NA,nrow=nrow(coefmin),ncol=nrow(predloc))
for(kk in 1:nrow(coefmin.sim)){
  coefmin.sim[kk,] = suppressWarnings(predict(Krig(lon.lat,coefmin[kk,]),predloc))
}
# maximum temperature
coefmax.sim <- matrix(NA,nrow=nrow(coefmax),ncol=nrow(predloc))
for(kk in 1:nrow(coefmax.sim)){
  coefmax.sim[kk,] = suppressWarnings(predict(Krig(lon.lat,coefmax[kk,]),predloc))
}


##########################################
##
## Simulations
##
##########################################
# simulation metadata
x_coords = unique(predloc.GK[,1]) * 1000 # convert back to meters (original form)
y_coords = unique(predloc.GK[,2]) * 1000 # convert back to meters (original form)
n.x_coords <- length(x_coords)
n.y_coords <- length(y_coords)
# elevation at each grid cell 
elev <- as.matrix(read.table("../1_setup/data/Salado-Abasin-elevation-meters.dat",header=F))
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
NT=2#00
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

# Gamma shape and scale
SH <- numeric(np)
SC <- array(data=NA,dim=c(nt.sim,np,NT))

for(k in 1:np) SH[k] <- gamma.shape(GAMMA[[k]])$alpha

for(i in 1:NT){
  for(k in 1:np){
    SC[,k,i] <- exp(apply(GAMMA[[k]]$coef*rbind(1,ct.sim,st.sim,ST1.sim[i,],ST2.sim[i,],ST3.sim[i,],ST4.sim[i,]),FUN=sum,MAR=2,na.rm=T))/SH[k]

  }
}

SH.sim <- suppressWarnings(predict(Krig(lon.lat,SH),predloc))


# define arrays for simulated series
SIMamt.sim <- SIMocc.sim <- SIMmax.sim <- SIMmin.sim <- array(dim=c(nt.sim,nrow(predloc),NT))

# TODO: remove this.
set.seed(1234)

# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
## occurrences
for(i in 1:NT){
  # initialize with climatology of Sept 30
  w2 <- suppressWarnings(grf(nrow(predloc),grid=predloc,cov.model="exponential",
            cov.pars=c(params[9,2],params[9,3]),nugget=params[9,1],mean=rep(0,nrow(predloc)),messages=FALSE))
  SIMocc.old <- (w2$data > 0) + 0

  tmin.points = apply(MN[da==30 & mo==9,],2,mean,na.rm=T)
  tmax.points = apply(MX[da==30 & mo==9,],2,mean,na.rm=T)
  SIMmin.old <- as.vector(predict(Krig(lon.lat,tmin.points),predloc))
  SIMmax.old <- as.vector(predict(Krig(lon.lat,tmax.points),predloc))

for(d in 1:nt.sim){
  # X.OCC <- cbind(POCC[,i],ct,st,ST1,ST2,ST3,ST4)
    require(microbenchmark)
    # microbenchmark(mu_occ=apply(rbind(1,SIMocc.old,ct.sim[d],st.sim[d],ST1.sim[i,d],ST2.sim[i,d],ST3.sim[i,d],ST4.sim[i,d])*coefocc.sim,2,sum,na.rm=T),
    #                w_occ=suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential", cov.pars=c(params[mo.sim[d],2],params[mo.sim[d],3]),nugget=params[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data))
  mu.occ <- apply(rbind(1,SIMocc.old,ct.sim[d],st.sim[d],ST1.sim[i,d],ST2.sim[i,d],ST3.sim[i,d],ST4.sim[i,d])*coefocc.sim,2,sum,na.rm=T)
  w.occ  <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params[mo.sim[d],2],params[mo.sim[d],3]),nugget=params[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
  SIMocc.sim[d,,i] <- ((mu.occ+w.occ) > 0) + 0
    # X.MIN <- X.MAX <- cbind(PMN[,i], PMX[,i], ct, st, OCC[,i], Rt, SMN1, SMN2, SMN3, SMN4, SMX1, SMX2, SMX3, SMX4)
  mu.min <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc.sim[d,,i],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmin.sim,2,sum,na.rm=T)
  mu.max <- apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc.sim[d,,i],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmax.sim,2,sum,na.rm=T)
  
  microbenchmark(mu_occ=apply(rbind(1,SIMocc.old,ct.sim[d],st.sim[d],ST1.sim[i,d],ST2.sim[i,d],ST3.sim[i,d],ST4.sim[i,d])*coefocc.sim,2,sum,na.rm=T),
                 mu_min=apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc.sim[d,,i],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmin.sim,2,sum,na.rm=T),
                 mu_max=apply(rbind(1,SIMmin.old,SIMmax.old,ct.sim[d],st.sim[d],SIMocc.sim[d,,i],Rt.sim[d],SMN1.sim[i,d],SMN2.sim[i,d],SMN3.sim[i,d],SMN4.sim[i,d],SMX1.sim[i,d],SMX2.sim[i,d],SMX3.sim[i,d],SMX4.sim[i,d])*coefmax.sim,2,sum,na.rm=T),
                 w_occ=suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential", cov.pars=c(params[mo.sim[d],2],params[mo.sim[d],3]),nugget=params[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data),
                 w_min=suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential", cov.pars=c(params.min[mo.sim[d],2],params.min[mo.sim[d],3]),nugget=params.min[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data),
                 w_max=suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential", cov.pars=c(params.max[mo.sim[d],2],params.max[mo.sim[d],3]),nugget=params.max[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data))
  
  w.min <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params.min[mo.sim[d],2],params.min[mo.sim[d],3]),nugget=params.min[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
  w.max <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params.max[mo.sim[d],2],params.max[mo.sim[d],3]),nugget=params.max[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
  
  SIMmin.sim[d,,i] <- signif(mu.min+w.min,digits=4)
  SIMmax.sim[d,,i] <- signif(mu.max+w.max,digits=4)

    while(min(SIMmax.sim[d,,i] - SIMmin.sim[d,,i]) < 0.5){
      w.min <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params.min[mo.sim[d],2],params.min[mo.sim[d],3]),nugget=params.min[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
      w.max <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params.max[mo.sim[d],2],params.max[mo.sim[d],3]),nugget=params.max[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
  
      SIMmin.sim[d,,i] <- signif(mu.min+w.min,digits=4)
      SIMmax.sim[d,,i] <- signif(mu.max+w.max,digits=4)
    }

  SIMocc.old <- SIMocc.sim[d,,i]
  SIMmax.old <- SIMmax.sim[d,,i]
  SIMmin.old <- SIMmin.sim[d,,i]

  ## amounts
  SC.sim <- suppressWarnings(predict(Krig(lon.lat,SC[d,,i]),predloc))
  
  w3 <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
            cov.pars=c(params[mo.sim[d],2],params[mo.sim[d],3]),nugget=params[mo.sim[d],1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
  SIMamt.sim[d,,i] <- signif(qgamma(pnorm(w3), shape=SH.sim, scale=SC.sim),digits=4)

  }
  setTxtProgressBar(pb, i)
}
  SIMamt.sim[SIMamt.sim<0.1]= 0
  SIMocc.sim[SIMamt.sim==0] = 0
  SIMamt.sim[SIMocc.sim==0] = 0

# estimate potential evapotranspiration
# set up for calculating reference evapotranspiration (ET)
  day.of.year <- yday(as.Date(paste(yr.sim,mo.sim,da.sim,sep="-")))
  q0 = array(data=NA,dim=dim(SIMamt.sim))
  predrad <- numeric(nrow(predloc))

  for(kk in 1:nrow(predloc)){ 
    predrad[kk] <- radians(predloc[kk,2])
    q0[,kk,1] <- extrat(i=day.of.year, lat=predrad[kk])$"ExtraTerrestrialSolarRadiationDaily"
  }
  if(NT > 1) {
    for(i in 2:NT) q0[,,i] <- q0[,,1]
  }

# set up for calculating reference evapotranspiration (ET)
  coef.a <- 0.001703
  coef.b <- 21.967919
  coef.c <- 0.083444
  coef.d <- 0.541066
# set up for calculating solar radiation (SRAD)
  bc.coef=c(A=0.69, B=0.02, C=2.12)

# difference between max and min temperatures
  dtr <- SIMmax.sim-SIMmin.sim
# we included a check in the wgen code to ensure max temp is always greater than min temp... so the next line is not necessary
  #dtr[dtr<0] <- NA
# average temperature
  tavg = (SIMmax.sim+SIMmin.sim)/2
# equation for reference ET
SIMet0.sim <- coef.a * 0.408 * q0 * (tavg + coef.b) * (dtr - (coef.c * SIMamt.sim)) ** coef.d
SIMet0.sim[is.na(SIMet0.sim)]=0 # R cant handle imaginary numbers (i.e. exponentiating negative bases)


# Estimating solar radiation.
# Computate Bristow-Campbell's delta temperature.
# do not consider tomorrow's minimum temperature 
# because this is hourly aggregate Met Service data
# dtr is dtemp
extraT <- array(suppressWarnings(extrat(i=matrix(dayOfYear(matrix(rep(sim.dates,nrow(predloc.GK)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F), lat=predrad)$ExtraTerrestrialSolarRadiationDaily),dim=dim(dtr))
SIMsrad.sim <- extraT * bc.coef[['A']] * (1 - exp(-bc.coef[['B']] * (dtr^bc.coef[['C']])))


# define dimensions for .nc file 
dimRNum <- ncdim_def(name="rnum", units="number", vals=RNum, unlim=TRUE, create_dimvar=TRUE, longname="Realization number")
dimX <- ncdim_def(name="x_coord", units="meters", vals=x_coords, unlim=FALSE, create_dimvar=TRUE, longname="longitude in Gauss-Krueger coordinates")
dimY <- ncdim_def(name="y_coord", units="meters", vals=y_coords, unlim=FALSE, create_dimvar=TRUE, longname="latitude in Gauss-Krueger coordinates")
dimTime <- ncdim_def(name="time", units="days since 1961-01-01", vals=Times, unlim = FALSE, create_dimvar=TRUE, calendar="standard", longname="Time in days since 1961-01-01")

tmax <- ncvar_def(name="tmax", units="degress_Celsius", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily maximum near-surface air temperature", prec="float", verbose=TRUE)
tmin <- ncvar_def(name="tmin", units="degress_Celsius", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily minimum near-surface air temperature", prec="float", verbose=TRUE)
prcp <- ncvar_def(name="prcp", units="mm day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily total rainfall", prec="float", verbose=TRUE)
srad <- ncvar_def(name="srad", units="Mjoules m-2 day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily solar radiation", prec="float", verbose=TRUE)
et0  <- ncvar_def(name="et0",  units="mm day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Reference evapotranspiration", prec="float", verbose=TRUE)
dec_lon <- ncvar_def(name="dec_lon", units="degrees_east", missval=mv, dim=list(dimY, dimX), longname="Longitude in decimal degrees", prec="float", verbose=TRUE)
dec_lat <- ncvar_def(name="dec_lat", units="degrees_north", missval=mv, dim=list(dimY, dimX), longname="Latitude in decimal degrees", prec="float", verbose=TRUE)
elevation <- ncvar_def(name="elevation", units="meters from sea level", missval=mv, dim=list(dimY, dimX), longname="Elevation in meters from sea level", prec="float", verbose=TRUE)

STs <- ncvar_def(name="ST", units="mm season-1", dim=list(dimRNum), missval=mv, longname="Seasonal total precipitation covariates", prec="float", verbose=TRUE)
SMXs <- ncvar_def(name="SMX", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average maximum temperature covariates", prec="float", verbose=TRUE)
SMNs <- ncvar_def(name="SMN", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average minimum temperature covariates", prec="float", verbose=TRUE)

# Create NetCDF file in current working directory
dd <- nc_create(filename="wgen-conditional-seasonal-gridded.nc",
  vars = list(tmax, tmin, prcp, srad, et0, dec_lon, dec_lat, elevation, STs, SMXs, SMNs),
  force_v4 = TRUE, verbose = TRUE)

# define arrays for simulated series
# array(dim=c(NT,nt.sim,n.y_coords,n.x_coords))
# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
for(i in 1:NT){
  for(d in 1:nt.sim){
    ncvar_put(dd, varid = "tmax", 
      vals = matrix(SIMmax.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
      start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
      verbose = TRUE)

    ncvar_put(dd, varid = "tmin", 
      vals = matrix(SIMmin.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
      start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
      verbose = TRUE)

    ncvar_put(dd, varid = "prcp", 
      vals = matrix(SIMamt.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
      start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
      verbose = TRUE)
      
    ncvar_put(dd, varid = "et0", 
      vals = matrix(SIMet0.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
      start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
      verbose = TRUE)

    ncvar_put(dd, varid = "srad",
      vals = matrix(SIMsrad.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
      start = c(i, d, 1,1 ), count = c(1, 1, n.y_coords, n.x_coords),
      verbose = TRUE)
  }
setTxtProgressBar(pb, i)
}

ncvar_put(dd, varid = "dec_lon",
  vals = matrix(predlocs[,"lon"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
  start = c(1, 1), count = c(n.y_coords, n.x_coords),
  verbose = TRUE)

ncvar_put(dd, varid = "dec_lat",
  vals = matrix(predlocs[,"lat"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
  start = c(1, 1), count = c(n.y_coords, n.x_coords),
  verbose = TRUE)

ncvar_put(dd, varid = "elevation",
  vals = as.matrix(elev),
  start = c(1, 1), count = c(n.y_coords, n.x_coords),
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
