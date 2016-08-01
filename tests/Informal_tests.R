rm(list=objects()); gc();

library(glmwgen)
require(sp)
projections <- list(
    # Gauss-Kruger
    'gk' = paste('+proj=tmerc +lat_0=-90 +lon_0=-60 +k=1 +x_0=5500000 +y_0=0',
                 '+ellps=intl +twogs84=-148,136,90,0,0,0,0 +units=m +no_defs'),
    # Latitude and longitude in decimal degrees.
    'latlon' = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
)

load(file='data/test/stations.RData')
load(file='data/test/climate.RData')

require(doMC)
# registerDoMC(detectCores())
registerDoMC(2)

stations_gk <- sp::spTransform(stations, CRS(projections$gk))

glmgen_fit <- calibrate.glmwgen(climate, stations)
glmgen_fit_gk <- calibrate.glmwgen(climate, stations_gk)

# unique_lon <- unique(coordinates(stations)[,1])
# unique_lat <- unique(coordinates(stations)[,2])
#
# unique_x <- unique(coordinates(stations_gk)[,1])
# unique_y <- unique(coordinates(stations_gk)[,2])
#
# grid <- expand.grid(x=seq(min(unique_lon), max(unique_lon), by = 0.5), y=seq(min(unique_lat), max(unique_lat), by = 0.5))
# grid_gk <- expand.grid(x=seq(min(unique_x), max(unique_x), by = 10000), y=seq(min(unique_y), max(unique_y), by = 10000))

grid <- read.csv('_workshop/1_setup/data/grids.txt', sep ='\t')
grid_gk <- SpatialPoints(grid[, c(1, 2)], proj4string = stations_gk@proj4string)
grid <- SpatialPoints(grid[, c(3, 4)], proj4string = stations@proj4string)

# plot(make_distance_grid(grid), pch = '.')

microbenchmark::microbenchmark({
    simulated_climate <- simulate(glmgen_fit, start_date = '2016-01-01', end_date = '2016-02-01')
}, times = 10)

# Unit: seconds
#      min       lq     mean   median       uq      max  neval
# 6.439614 6.648005  6.72933 6.730807 6.769429 7.017718     10

# microbenchmark::microbenchmark({
#     simulated_climate <- simulate(glmgen_fit, start_date = '2016-01-01', end_date = '2016-02-01', control = glmwgenSimulationControl(grf_method = 'svd'))
# }, times = 10)


system.time(simulated_climate_latlon <- simulate(glmgen_fit, start_date = '2016-06-01', end_date = '2016-06-15',
                                          simulation_locations = grid, seed = 123))
system.time(simulated_climate_gk <- simulate(glmgen_fit_gk, start_date = '2016-06-01', end_date = '2016-06-15',
                                          simulation_locations = grid_gk, seed = 123))
# system.time(simulated_climate_gk2 <- simulate(glmgen_fit_gk, start_date = '2016-01-01', end_date = '2016-02-01',
#                                           simulation_locations = grid_gk, seed = 1234))

# system.time(simulated_climate_chol <- simulate(glmgen_fit, start_date = '2016-01-01', end_date = '2016-06-01',
#                                                simulation_locations = grid,
#                                                control = glmwgenSimulationControl(random_fields_method = 'chol')))

# fields::quilt.plot(coordinates(grid_gk), simulated_climate_gk[[1]]$prcp[1,])
#
# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$tx[12,])
#
# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$prcp[2,])
# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$tx[2,])
# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$prcp[12,])
# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$tx[12,])

# fields::quilt.plot(grid, simulated_climate_latlon[[1]]$tx[14,])
# fields::quilt.plot(grid_gk, simulated_climate_gk2[[1]]$tx[1,])


for(i in 1:15) {
    fields::quilt.plot(coordinates(grid_gk), simulated_climate_gk[[1]]$prcp[i,], main = paste('GK Day', i))
    fields::quilt.plot(coordinates(grid), simulated_climate_latlon[[1]]$prcp[i,], main = paste('LatLong Day', i))
}


fields::quilt.plot(grid_gk, simulated_climate2[[1]]$tx[2,])
fields::quilt.plot(grid_gk, simulated_climate[[1]]$tx[2,])

fields::quilt.plot(grid_gk, simulated_climate[[1]]$tn[1,])
fields::quilt.plot(grid_gk, simulated_climate[[1]]$prcp[1,])

fields::quilt.plot(grid_gk, simulated_climate[[1]]$tx[15,])
fields::quilt.plot(grid_gk, simulated_climate[[1]]$tn[15,])
fields::quilt.plot(grid_gk, simulated_climate[[1]]$prcp[15,])
# fields::quilt.plot(grid, simulated_climate[[1]]$prcp[1,])
# fields::quilt.plot(grid, simulated_climate_chol[[1]]$tx[1,])
