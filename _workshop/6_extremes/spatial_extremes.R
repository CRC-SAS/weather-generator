  # --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://lib.stat.cmu.edu/R/CRAN/"), error=traceback)

# only need ismev package for computing GEV parameters
# and fields package for plotting the spatial maps
# and rgdal to read in the shapefile for A basin
list.of.packages <- c("ismev", "fields", "rgdal")  

for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack)
    require(pack, character.only = TRUE)
  }
}

rm(list.of.packages, pack)

## set the working directory -- this is wherever you put the folder from the zip file
setwd("/Users/andrew/Desktop/workshop-BA/6_extremes")
## load the setup data
load("../1_setup/wgen-setup.RData")

predloc <- read.table("../1_setup/data/grids.txt",header=T)[,3:4]
predloc.GK <- read.table("../1_setup/data/grids.txt",header=T)[,1:2]
cuenca_A_ll <- readOGR("../1_setup/data/cuenca", layer="Salado_A_latlon")
cuenca_A_gk <- readOGR("../1_setup/data/cuenca", layer="Salado_A_GK")

elev <- stns[,"elev"]
elev.predloc <- as.matrix(read.table("../1_setup/data/Salado-Abasin-elevation-meters.dat",header=F))

sea.start1 <- which(mo==1 & da==1)
sea.end1  <- which(mo==3 & da==31)
sea.start2 <- which(mo==4 & da==1)
sea.end2  <- which(mo==6 & da==30)
sea.start3 <- which(mo==7 & da==1)
sea.end3  <- which(mo==9 & da==30)
sea.start4 <- which(mo==10 & da==1)
sea.end4  <- which(mo==12 & da==31)

PP.ex1 <- PP.ex2 <- PP.ex3 <- PP.ex4 <- PP.na1 <- PP.na2 <- PP.na3 <- PP.na4 <- 
MX.ex1 <- MX.ex2 <- MX.ex3 <- MX.ex4 <- MX.na1 <- MX.na2 <- MX.na3 <- MX.na4 <-
MN.ex1 <- MN.ex2 <- MN.ex3 <- MN.ex4 <- MN.na1 <- MN.na2 <- MN.na3 <- MN.na4 <- 
                                        matrix(NA,nrow=nyr,ncol=np)

for(k in 1:np){
  for(j in 1:nyr){
    
    # these if statements ensure the values returned are not Inf or -Inf
    
    if(sum(is.na(PP[sea.start1[j]:sea.end1[j],k])) == length(PP[sea.start1[j]:sea.end1[j],k])){
      next
    }else{
      PP.ex1[j,k] <- max(PP[sea.start1[j]:sea.end1[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      PP.na1[j,k] <- sum(is.na(PP[sea.start1[j]:sea.end1[j],k]))
    }
    
    if(sum(is.na(PP[sea.start2[j]:sea.end2[j],k])) == length(PP[sea.start2[j]:sea.end2[j],k])){
      next
    }else{
      PP.ex2[j,k] <- max(PP[sea.start2[j]:sea.end2[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      PP.na2[j,k] <- sum(is.na(PP[sea.start2[j]:sea.end2[j],k]))
    }
    
    if(sum(is.na(PP[sea.start3[j]:sea.end3[j],k])) == length(PP[sea.start3[j]:sea.end3[j],k])){
      next 
    }else{
      PP.ex3[j,k] <- max(PP[sea.start3[j]:sea.end3[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      PP.na3[j,k] <- sum(is.na(PP[sea.start3[j]:sea.end3[j],k]))
    }
    
    if(sum(is.na(PP[sea.start4[j]:sea.end4[j],k])) == length(PP[sea.start4[j]:sea.end4[j],k])){
      next 
    }else{
      PP.ex4[j,k] <- max(PP[sea.start4[j]:sea.end4[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      PP.na4[j,k] <- sum(is.na(PP[sea.start4[j]:sea.end4[j],k]))
    }
    
    if(sum(is.na(MX[sea.start1[j]:sea.end1[j],k])) == length(MX[sea.start1[j]:sea.end1[j],k])){
      next 
    }else{
      MX.ex1[j,k] <- max(MX[sea.start1[j]:sea.end1[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MX.na1[j,k] <- sum(is.na(MX[sea.start1[j]:sea.end1[j],k]))
    }

    if(sum(is.na(MX[sea.start2[j]:sea.end2[j],k])) == length(MX[sea.start2[j]:sea.end2[j],k])){
      next 
    }else{
      MX.ex2[j,k] <- max(MX[sea.start2[j]:sea.end2[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MX.na2[j,k] <- sum(is.na(MX[sea.start2[j]:sea.end2[j],k]))
    }

    if(sum(is.na(MX[sea.start3[j]:sea.end3[j],k])) == length(MX[sea.start3[j]:sea.end3[j],k])){
      next
    }else{
      MX.ex3[j,k] <- max(MX[sea.start3[j]:sea.end3[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MX.na3[j,k] <- sum(is.na(MX[sea.start3[j]:sea.end3[j],k]))
    }

    if(sum(is.na(MX[sea.start4[j]:sea.end4[j],k])) == length(MX[sea.start4[j]:sea.end4[j],k])){
      next 
    }else{
      MX.ex4[j,k] <- max(MX[sea.start4[j]:sea.end4[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MX.na4[j,k] <- sum(is.na(MX[sea.start4[j]:sea.end4[j],k]))
    }

    if(sum(is.na(MN[sea.start1[j]:sea.end1[j],k])) == length(MN[sea.start1[j]:sea.end1[j],k])){
      next 
    }else{
      MN.ex1[j,k] <- min(MN[sea.start1[j]:sea.end1[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MN.na1[j,k] <- sum(is.na(MN[sea.start1[j]:sea.end1[j],k]))
    }

    if(sum(is.na(MN[sea.start2[j]:sea.end2[j],k])) == length(MN[sea.start2[j]:sea.end2[j],k])){
      next 
    }else{
      MN.ex2[j,k] <- min(MN[sea.start2[j]:sea.end2[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MN.na2[j,k] <- sum(is.na(MN[sea.start2[j]:sea.end2[j],k]))
    }

    if(sum(is.na(MN[sea.start3[j]:sea.end3[j],k])) == length(MN[sea.start3[j]:sea.end3[j],k])){
      next 
    }else{
      MN.ex3[j,k] <- min(MN[sea.start3[j]:sea.end3[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MN.na3[j,k] <- sum(is.na(MN[sea.start3[j]:sea.end3[j],k]))
    }

    if(sum(is.na(MN[sea.start4[j]:sea.end4[j],k])) == length(MN[sea.start4[j]:sea.end4[j],k])){
      next 
    }else{
      MN.ex4[j,k] <- min(MN[sea.start4[j]:sea.end4[j],k],na.rm=T)
      # also want to remove seasons where there is excessive missing data
      MN.na4[j,k] <- sum(is.na(MN[sea.start4[j]:sea.end4[j],k]))
    }

  }
}

# mask the seasons with excessive missing data as NA 
PP.ex1[PP.na1 > 60] = NA
PP.ex2[PP.na2 > 60] = NA 
PP.ex3[PP.na3 > 60] = NA 
PP.ex4[PP.na4 > 60] = NA 

MX.ex1[MX.na1 > 60] = NA
MX.ex2[MX.na2 > 60] = NA 
MX.ex3[MX.na3 > 60] = NA 
MX.ex4[MX.na4 > 60] = NA 

MN.ex1[MN.na1 > 60] = NA 
MN.ex2[MN.na2 > 60] = NA 
MN.ex3[MN.na3 > 60] = NA 
MN.ex4[MN.na4 > 60] = NA 

# we will now estimate return levels for 100 year return period
return.period <- 100

# we will now produce spatial maps of JFM extreme precip
# define vectors for GEV location, shape, and scale parameters 
# these are spatially varying, so we can produce spatial map of
# 100-year (or any other return period) extreme precip...
# this is a stationary approach, which does not consider climate change
locs.JFM <- numeric(np)
shps.JFM <- numeric(np)
scls.JFM <- numeric(np)
for(k in 1:np){
locs.JFM[k] <- gev.fit(PP.ex1[,k][!is.na(PP.ex1[,k])])$mle[1]
scls.JFM[k] <- gev.fit(PP.ex1[,k][!is.na(PP.ex1[,k])])$mle[2]
shps.JFM[k] <- gev.fit(PP.ex1[,k][!is.na(PP.ex1[,k])])$mle[3]
}
# estimate the parameters on the grid
locmat.JFM <- predict(Krig(x=lon.lat,Y=locs.JFM,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
sclmat.JFM <- predict(Krig(x=lon.lat,Y=scls.JFM,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
shpmat.JFM <- predict(Krig(x=lon.lat,Y=shps.JFM,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
# calculate using the GEV formula
return.level.JFM <- locmat.JFM + sclmat.JFM*((-log10(1-1/return.period))^-shpmat.JFM -1)/shpmat.JFM

# we will now produce spatial maps of AMJ extreme precip
# define vectors for GEV location, shape, and scale parameters 
# these are spatially varying, so we can produce spatial map of
# 100-year (or any other return period) extreme precip...
# this is a stationary approach, which does not consider climate change
locs.AMJ <- numeric(np)
shps.AMJ <- numeric(np)
scls.AMJ <- numeric(np)
for(k in 1:np){
locs.AMJ[k] <- gev.fit(PP.ex2[,k][!is.na(PP.ex2[,k])])$mle[1]
scls.AMJ[k] <- gev.fit(PP.ex2[,k][!is.na(PP.ex2[,k])])$mle[2]
shps.AMJ[k] <- gev.fit(PP.ex2[,k][!is.na(PP.ex2[,k])])$mle[3]
}
# estimate the paramters on the grid
locmat.AMJ <- predict(Krig(x=lon.lat,Y=locs.AMJ,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
sclmat.AMJ <- predict(Krig(x=lon.lat,Y=scls.AMJ,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
shpmat.AMJ <- predict(Krig(x=lon.lat,Y=shps.AMJ,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
# calculate using the GEV formula
return.level.AMJ <- locmat.AMJ + sclmat.AMJ*((-log10(1-1/return.period))^-shpmat.AMJ -1)/shpmat.AMJ

# we will now produce spatial maps of JAS extreme precip
# define vectors for GEV location, shape, and scale parameters 
# these are spatially varying, so we can produce spatial map of
# 100-year (or any other return period) extreme precip...
# this is a stationary approach, which does not consider climate change
locs.JAS <- numeric(np)
shps.JAS <- numeric(np)
scls.JAS <- numeric(np)
for(k in 1:np){
locs.JAS[k] <- gev.fit(PP.ex3[,k][!is.na(PP.ex3[,k])])$mle[1]
scls.JAS[k] <- gev.fit(PP.ex3[,k][!is.na(PP.ex3[,k])])$mle[2]
shps.JAS[k] <- gev.fit(PP.ex3[,k][!is.na(PP.ex3[,k])])$mle[3]
}
# estimate the parameters on the grid
locmat.JAS <- predict(Krig(x=lon.lat,Y=locs.JAS,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
sclmat.JAS <- predict(Krig(x=lon.lat,Y=scls.JAS,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
shpmat.JAS <- predict(Krig(x=lon.lat,Y=shps.JAS,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
# calculate using the GEV formula
return.level.JAS <- locmat.JAS + sclmat.JAS*((-log10(1-1/return.period))^-shpmat.JAS -1)/shpmat.JAS

# we will now produce spatial maps of OND extreme precip
# define vectors for GEV location, shape, and scale parameters 
# these are spatially varying, so we can produce spatial map of
# 100-year (or any other return period) extreme precip...
# this is a stationary approach, which does not consider climate change
locs.OND <- numeric(np)
shps.OND <- numeric(np)
scls.OND <- numeric(np)
for(k in 1:np){
locs.OND[k] <- gev.fit(PP.ex4[,k][!is.na(PP.ex4[,k])])$mle[1]
scls.OND[k] <- gev.fit(PP.ex4[,k][!is.na(PP.ex4[,k])])$mle[2]
shps.OND[k] <- gev.fit(PP.ex4[,k][!is.na(PP.ex4[,k])])$mle[3]
}
# estimate the parameters on the grid
locmat.OND <- predict(Krig(x=lon.lat,Y=locs.OND,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
sclmat.OND <- predict(Krig(x=lon.lat,Y=scls.OND,Z=elev,m=1),predloc, Z=as.vector(t(elev.predloc)))
shpmat.OND <- predict(Krig(x=lon.lat,Y=shps.OND,Z=elev, m=1),predloc, Z=as.vector(t(elev.predloc)))
# calculate using the GEV formula
return.level.OND <- locmat.OND + sclmat.OND*((-log10(1-1/return.period))^-shpmat.OND -1)/shpmat.OND




par(mfrow=c(2,2), mar=c(2,2,2,2))
# JFM
quilt.plot(predloc.GK, return.level.JFM, main=paste(return.period,"Year Return Level for JFM"),axes=F)
box()
plot(cuenca_A_gk, add=T, lwd=2)
# AMJ
quilt.plot(predloc.GK, return.level.AMJ, main=paste(return.period,"Year Return Level for AMJ"),axes=F)
box()
plot(cuenca_A_gk, add=T, lwd=2)
# JAS
quilt.plot(predloc.GK, return.level.JAS, main=paste(return.period,"Year Return Level for JAS"),axes=F)
box()
plot(cuenca_A_gk, add=T, lwd=2)
# OND
quilt.plot(predloc.GK, return.level.OND, main=paste(return.period,"Year Return Level for OND"),axes=F)
box()
plot(cuenca_A_gk, add=T, lwd=2)


