  # --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://lib.stat.cmu.edu/R/CRAN/"), error=traceback)

# only need ismev package for computing GEV parameters
# and fields package for plotting the spatial maps
# and rgdal to read in the shapefile for A basin
list.of.packages <- c("geoR", "ncdf4", "fields", "lubridate", "MASS")  

for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack)
    require(pack, character.only = TRUE)
  }
}

rm(list.of.packages, pack)

setwd("/Users/andrew/Desktop/workshop-BA/5_metamodel")
PP <- read.table("../1_setup/data/Salado_prcp.dat",header=T)
yr <- year(PP[,1])
mo <- month(PP[,1])
da <- day(PP[,1])
uyr <- unique(yr)
nyr <- length(uyr)
rm(PP)

dd.in <- nc_open("../1_setup/data/2015-11-11-BigPampa-observed-interpolation-1961-2013.nc")
dd.out <- nc_open("../1_setup/data/mike_output.nc")
dd.et <- nc_open("../1_setup/data/mike_output_et.nc")

# setwd("../")
prcp <- ncvar_get(dd.in, "prcp")
tmax <- ncvar_get(dd.in, "tmax")
tmin <- ncvar_get(dd.in, "tmin")
srad <- ncvar_get(dd.in, "srad")
et0  <- ncvar_get(dd.in, "et0")

napa <- ncvar_get(dd.out, "wdepth")
elev <- ncvar_get(dd.out, "elevation")

et <- ncvar_get(dd.et, "wdepth")

lon.GK <- ncvar_get(dd.out, "x_coord")/1000
lat.GK <- ncvar_get(dd.out, "y_coord")/1000

nc_close(dd.in)
nc_close(dd.out)
nc_close(dd.et)

predloc <- read.table("../1_setup/data/grids.txt",header=T)
lon <- matrix(predloc[,3],nrow=dim(napa)[2],ncol=dim(napa)[3],byrow=T)
lat <- matrix(predloc[,4],nrow=dim(napa)[2],ncol=dim(napa)[3],byrow=T)

# number of days in data set
tt <- nrow(napa)

# 12 months per year, nyr years in data set
tt.month <- 12*nyr
yr.month <- rep(uyr,each=12)
mo.month <- rep(1:12,nyr)

nx <- dim(napa)[2]
ny <- dim(napa)[3]

xx.id <- seq(from=0,to=nx,by=2)
yy.id <- seq(from=0,to=ny,by=2)

boundary.mask <- basin.mask <- matrix(data=NA,nrow=dim(napa)[2],ncol=dim(napa)[3])

for(xx in 1:nx){
	for(yy in 1:ny){
		if(is.na(mean(napa[2:tt,xx,yy])) | (mean(napa[2:tt,xx,yy]) == -1)){
			if(!is.na(mean(napa[2:tt,xx,yy])) & mean(napa[2:tt,xx,yy]) == -1){
				boundary.mask[xx,yy] = -1
				basin.mask[xx,yy] = 1
			}else{
			next
			}
		}
		basin.mask[xx,yy] = 1
	}
}

new.napa <- array(data=NA,dim=dim(napa)[c(1:3)])

for(k in 1:tt){
	print(k)
	look <- image.smooth(napa[k,,], dx=5, dy=5, theta=5)
	new.napa[k,,] <- look$z 
	new.napa[k,,][is.na(basin.mask)] = NA 
	new.napa[k,,][boundary.mask==-1] = NA
}

napa.m.s <- array(data=NA,dim=c(tt.month,dim(new.napa)[2:3]))

for(k in 1:tt.month){
	print(k)
	month.id <- yr==yr.month[k] & mo==mo.month[k]

	napa.m.s[k,,] <- apply(new.napa[month.id,,],2:3,mean)
}

prcp.m <- tmax.m <- tmin.m <- srad.m <- et0.m <- et.m <- napa.m <- 
			array(data=NA,dim=c(tt.month,dim(napa)[2:3]))

for(k in 1:tt.month){
	print(k)
	month.id <- yr==yr.month[k] & mo==mo.month[k]
	prcp.m[k,,] <- apply(prcp[month.id,,],c(2:3),mean)
	tmax.m[k,,] <- apply(tmax[month.id,,],c(2:3),mean)
	tmin.m[k,,] <- apply(tmin[month.id,,],c(2:3),mean)
	srad.m[k,,] <- apply(srad[month.id,,],c(2:3),mean)
	et0.m[k,,]  <- apply(et0[month.id,,],c(2:3),mean)
	et.m[k,,]   <- apply(et[month.id,,],c(2:3),mean)

	napa.m[k,,] <- apply(napa[month.id,,],c(2:3),mean)
}

girasol.y <- readRDS("../1_setup/data/landUse-girasol.rds")
maize.y <- readRDS("../1_setup/data/landUse-maize.rds")
pasturas.y <- readRDS("../1_setup/data/landUse-pasturas.rds")
soja1.y <- readRDS("../1_setup/data/landUse-soja1.rds")
trigo.y <- readRDS("../1_setup/data/landUse-trigo.rds")
trigosoja2.y <- readRDS("../1_setup/data/landUse-trigosoja2.rds")

girasol.m <- maize.m <- pasturas.m <- soja1.m <- trigo.m <- trigosoja2.m <- array(data=NA,dim=c(tt.month,dim(napa)[2:3]))

for(k in 1:tt.month){
	print(k)
	index.id <- (yr.month[k] - 1960)
	girasol.m[k,,]    <- t(girasol.y[index.id,,])
	maize.m[k,,]      <- t(maize.y[index.id,,])
	pasturas.m[k,,]   <- t(pasturas.y[index.id,,])
	soja1.m[k,,]      <- t(soja1.y[index.id,,])
	trigo.m[k,,]      <- t(trigo.y[index.id,,])
	trigosoja2.m[k,,] <- t(trigosoja2.y[index.id,,])
}

save.image("metamodel-setup.RData")
