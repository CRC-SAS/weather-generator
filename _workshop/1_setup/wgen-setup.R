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
setwd("/home/federico/Desarrollo/RProjects/Wth Simulation/1_setup/")

## load the setup data
# station metadata
stns = read.table('data/Salado_metadata.dat',header=T)
## Lon.Lat && distance matrix
lon.lat = stns[,3:4]; colnames(lon.lat)=c("lon","lat")
dist.mat <- rdist.earth(lon.lat,miles=F)
diag(dist.mat) <- 0
# read in data
PP2 = read.table('data/Salado_prcp.dat',header=T)
MX2 = read.table('data/Salado_tmax.dat',header=T)
MN2 = read.table('data/Salado_tmin.dat',header=T)
# redefine data sets without dates for processing
PP = PP2[,2:ncol(PP2)]
MX = MX2[,2:ncol(MX2)]
MN = MN2[,2:ncol(MN2)]
# define number of days, number of stations
nt = nrow(PP)
np = ncol(PP)

## year, month, day vectors
yr <- year(PP2[,1])
mo <- month(PP2[,1]) 
da <- day(PP2[,1])
## years of record
uyr <- unique(yr)
nyr <- length(uyr)
firstyear <- uyr[1]
lastyear  <- uyr[length(uyr)]

## define LEAP YEARS
leapID = which(trunc(uyr/4)==uyr/4)
yearID = rep(365,nyr); yearID[leapID] = 366

##################################################
## Covariates
##################################################
# cos(t) and sin(t)
ct <- st <- c()
for(i in 1:nyr){
  ct <- c(ct,cos((2*pi*(1:yearID[i]))/yearID[i]))
  st <- c(st,sin((2*pi*(1:yearID[i]))/yearID[i]))
}

# OCC  = occurrence (where recorded precip is at least equal to 0.1mm)
OCC <- (PP >= 0.1) + 0
# POCC = previous day's occurrence
POCC <- OCC 
POCC[1,] <- NA
POCC[2:nrow(POCC),] <- OCC[1:(nrow(OCC)-1),]

# PPI holds only positive precipitation values 
PPI <- PP # PPI is precip intensity, so cut out zeros
PPI[PPI < 0.1] <- NA
# PMN = previous day's minimum temperature
PMN <- MN
PMN[1,] = NA
PMN[2:nrow(PMN),] <- MN[1:(nrow(MN)-1),]
# PMX = previous day's maximum temperature
PMX <- MX
PMX[1,] = NA
PMX[2:nrow(PMX),] <- MX[1:(nrow(MX)-1),]
# Rt = linear trend from -1 to +1 to account for temperature trends
Rt = seq(from = -1, to = 1, length.out = nrow(OCC))

# number of days
nt <- nrow(OCC)

save.image("wgen-setup.RData")
