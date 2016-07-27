rm(list=ls()); invisible(gc());
# --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://mirror.fcaglp.unlp.edu.ar/CRAN/"), error=traceback)

list.of.packages <- c("ncdf4","MASS","fields","geoR","lubridate","sirad","dplyr", "foreach", "abind")

for (pack in list.of.packages) {
    if (!require(pack, character.only = TRUE)) {
        install.packages(pack)
        require(pack, character.only = TRUE)
    }
}

rm(list.of.packages, pack)

source('R/optimize.R')
source('R/fit.R')
source('R/simulate.R')
source('R/seasonal.R')

# Stations climate.
# climate <- read.table('1_setup/data/Salado.dat')
climate <- read.table('data/Pampas.dat')
# climate$date <- as.Date(climate$date)
climate$date <- strptime(climate$date, format='%d/%m/%Y')
climate$date <- as.Date(climate$date)

# station metadata
# stns = read.table('1_setup/data/Salado_metadata.dat',header=T) %>% arrange(omm_id)
stns <- read.table('data/Pampas_metadata.dat', header=T) %>% arrange(omm_id)
## Lon.Lat && distance matrix
# lon.lat = stns[,3:4];
# 
# colnames(lon.lat) <- c("lon","lat")
# rownames(lon.lat) <- stns$omm_id
stns$id <- stns$omm_id

stations <- SpatialPointsDataFrame(coords=stns[, c('lon_dec', 'lat_dec')],
                                   data=stns,
                                   proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

# save(stations, file='data/stations.RData')
# save(climate, file='data/climate.RData')
# predlocs <- read.table("1_setup/data/grids.txt",header=T)
# predloc.GK <- predlocs[,1:2]/1000 # convert from meters to KILOMETERS
# predloc <- predlocs[,3:4]
# rm(predlocs)


if(require(doMC)) {
    registerDoMC(3)
}

m <- calibrate.glmgen(climate, stations)

# for(i in c(-1, 0, 1, 2, 5, 10)) {
#     Rt <- i
#     sim <- simulate.glmgen(m, start_date='2016-01-01', end_date='2025-12-31',
#                            simulation_locations = coordinates(stations[stations$id %in% c(87374, 87544, 87688),]),
#                            control = glmgenSimulationControl(Rt = Rt))
#     
#     plot(sim[[1]]$tn[, 1], type='l', main = paste0('Mean: ', mean(sim[[1]]$tn[, 1]), ', Rt: ', Rt))
#     abline(lm(sim[[1]]$tn[, 1] ~ seq_along(sim[[1]]$tn[, 1])), col='red', lwd=2)
# }

# t <- foreach(i=1:5, .combine=partial(abind, along=0.5), .multicombine=T) %do% {
#     return(matrix(array(dim=c(3, 10)), array(dim=c(3, 10)), ncol=2, dimnames=list(c(), c('tn', 'tx'))))
# }


gen_climate <- simulate.glmgen(m, start_date='2016-01-01', end_date='2025-12-31',
                               simulation_locations = coordinates(stations[stations$id %in% c(87374, 87544, 87688),]),
                               control = glmgenSimulationControl(Rt = 2),
                               n_realizations = 10)

# climate_table <- foreach(n_realization=1:length(gen_climate), .combine=rbind, .multicombine=T) %dopar% {
#     realization_table <- c()
#     for(v_name in names(gen_climate[[n_realization]])) {
#         d <- dim(gen_climate[[n_realization]][[v_name]])
#         var_climate <- matrix(ncol = 3)
#         for(wth_station_idx in 1:d[2]) {
#             wth_station <- c(87374, 87544, 87688)[wth_station_idx]
#             
#             var_climate <- rbind(var_climate, cbind(seq(from=1, to=d[1]),
#                                                     wth_station,
#                                                     gen_climate[[n_realization]][[v_name]][, wth_station_idx]))
#         }
#         var_climate <- na.omit(var_climate)
#         colnames(var_climate) <- c('day', 'station', v_name)
#         # station_climate <- matrix(station_climate, nrow = d[1], ncol = length(gen_climate[[n_realization]]), byrow = F)
#         realization_table <- cbind(realization_table, var_climate)
#     }
#     return(cbind(n_realization, realization_table[, unique(colnames(realization_table))]))
#     # return(matrix(realization_table, ))
# }

save(climate_table, file='climate_table.RData')

# paste(sapply(gen_climate, FUN=function(x) round(mean(x$tn[, 1], na.rm = T), 2)), collapse = ', ')
# 
# plot(gen_climate[[1]]$tn[, 1], type='l', main = paste0('Mean: ', mean(gen_climate[[1]]$tn[, 1])))
# abline(lm(gen_climate[[1]]$tn[, 1] ~ seq_along(gen_climate[[1]]$tn[, 1])), col='red', lwd=2)

# save(sim, file='data/ref_Guille.RData')


# rm(models); invisible(gc());

simulation_dates <- seq.Date(from=as.Date('2016-01-01'), to=as.Date('2025-12-31'), by='day')
simulation_locations <- coordinates(stations[stations$id %in% c(87374, 87544, 87688),])
NT <- 10




# estimate potential evapotranspiration
# set up for calculating reference evapotranspiration (ET)
day_of_year <- yday(simulation_dates)
q0 <- array(data=NA, dim=dim(gen_climate[[1]]$prcp))
predrad <- numeric(nrow(simulation_locations))

for(kk in 1:nrow(simulation_locations)){ 
    predrad[kk] <- radians(simulation_locations[kk,2])
    q0[,kk] <- extrat(i=day_of_year, lat=predrad[kk])$"ExtraTerrestrialSolarRadiationDaily"
}

# set up for calculating reference evapotranspiration (ET)
coef.a <- 0.001703
coef.b <- 21.967919
coef.c <- 0.083444
coef.d <- 0.541066
# set up for calculating solar radiation (SRAD)
bc.coef=c(A=0.69, B=0.02, C=2.12)

for(i in 1:NT) {
    # difference between max and min temperatures
    dtr <- gen_climate[[i]]$tx - gen_climate[[i]]$tn
    # we included a check in the wgen code to ensure max temp is always greater than min temp... so the next line is not necessary
    #dtr[dtr<0] <- NA
    # average temperature
    tavg <-  (gen_climate[[i]]$tx + gen_climate[[i]]$tn) / 2
    # equation for reference ET
    # gen_climate[[i]]$et0 <- coef.a * 0.408 * q0 * (tavg + coef.b) * (dtr - (coef.c * gen_climate[[i]]$prcp)) ** coef.d
    # gen_climate[[i]]$et0[is.na(gen_climate[[i]]$et0)] <- 0 # R cant handle imaginary numbers (i.e. exponentiating negative bases)
    
    # Estimating solar radiation.
    # Computate Bristow-Campbell's delta temperature.
    # do not consider tomorrow's minimum temperature 
    # because this is hourly aggregate Met Service data
    # dtr is dtemp
    # extraT <- array(suppressWarnings(extrat(i=matrix(dayOfYear(matrix(rep(sim.dates,nrow(predloc.GK)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F), lat=predrad)$ExtraTerrestrialSolarRadiationDaily),dim=dim(dtr))
    gen_climate[[i]]$srad <- q0 * bc.coef[['A']] * (1 - exp(-bc.coef[['B']] * (dtr^bc.coef[['C']])))
}

# quilt.plot(predloc.GK, (gen_climate[[1]]$tx - gen_climate[[1]]$tn)[2, ])
quilt.plot(predloc.GK, gen_climate[[1]]$srad[1, ], main='Srad')
quilt.plot(predloc.GK, gen_climate[[1]]$et0[1, ], main='Et0')
quilt.plot(predloc.GK, gen_climate[[1]]$prcp[1, ], main='Prcp')

quilt.plot(predloc.GK, gen_climate[[1]]$srad[35, ], main='Srad')
# quilt.plot(predloc.GK, gen_climate[[1]]$et0[30, ], main='Et0')
quilt.plot(predloc.GK, gen_climate[[1]]$prcp[35, ], main='Prcp')


# define dimensions for .nc file 
dimRNum <- ncdim_def(name="rnum", units="number", vals=RNum, unlim=TRUE, create_dimvar=TRUE, longname="Realization number")
dimX <- ncdim_def(name="x_coord", units="meters", vals=x_coords, unlim=FALSE, create_dimvar=TRUE, longname="longitude in Gauss-Krueger coordinates")
dimY <- ncdim_def(name="y_coord", units="meters", vals=y_coords, unlim=FALSE, create_dimvar=TRUE, longname="latitude in Gauss-Krueger coordinates")
dimTime <- ncdim_def(name="time", units="days since 1961-01-01", vals=julian(simulation_dates$date, origin = as.Date("1961-01-01")), unlim = FALSE, create_dimvar=TRUE, calendar="standard", longname="Time in days since 1961-01-01")

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
    for(d in 1:nrow(simulation_dates)){
        ncvar_put(dd, varid = "tmax", 
                  vals = matrix(gen_climate[[i]]$tx[d,], nrow=n.y_coords, ncol=n.x_coords, byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords))
        
        ncvar_put(dd, varid = "tmin", 
                  vals = matrix(gen_climate[[i]]$tn[d,], nrow=n.y_coords, ncol=n.x_coords, byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords))
        
        ncvar_put(dd, varid = "prcp", 
                  vals = matrix(gen_climate[[i]]$prcp[d,], nrow=n.y_coords, ncol=n.x_coords, byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords))
        
        ncvar_put(dd, varid = "et0", 
                  vals = matrix(gen_climate[[i]]$et0[d,], nrow=n.y_coords, ncol=n.x_coords, byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords))
        
        ncvar_put(dd, varid = "srad",
                  vals = matrix(gen_climate[[i]]$srad[d,], nrow=n.y_coords, ncol=n.x_coords, byrow=T),
                  start = c(i, d, 1,1 ), count = c(1, 1, n.y_coords, n.x_coords))
    }
    setTxtProgressBar(pb, i)
}

ncvar_put(dd, varid = "dec_lon",
          vals = matrix(predloc[,"lon"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
          start = c(1, 1), count = c(n.y_coords, n.x_coords))

ncvar_put(dd, varid = "dec_lat",
          vals = matrix(predloc[,"lat"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
          start = c(1, 1), count = c(n.y_coords, n.x_coords))

# ncvar_put(dd, varid = "elevation",
#           vals = as.matrix(elev),
#           start = c(1, 1), count = c(n.y_coords, n.x_coords),
#           verbose = TRUE)

# ncvar_put(dd, varid = "ST",
#           vals = apply(ST4.sim, 1, unique),
#           start = 1, count = NT,
#           verbose = TRUE)
# 
# ncvar_put(dd, varid = "SMX", 
#           vals = apply(SMX4.sim, 1, unique),
#           start = 1, count = NT,
#           verbose = TRUE)
# 
# ncvar_put(dd, varid = "SMN",
#           vals = apply(SMN4.sim, 1, unique),
#           start = 1, count = NT, 
#           verbose = TRUE)

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
