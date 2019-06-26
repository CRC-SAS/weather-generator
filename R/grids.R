make_distance_grid <- function(fit_stations, simulation_locations) {
    is_projected <- sp::is.projected(fit_stations)

    simulation_locations <- sp::coordinates(simulation_locations)
    stations_locations <- sp::coordinates(fit_stations)
    bounds <- sp::bbox(rbind(simulation_locations, stations_locations))

    simulation_grid <- matrix(nrow = nrow(simulation_locations), ncol = 2)
    stations_grid <- matrix(nrow = nrow(stations_locations), ncol = 2)

    for(point_idx in 1:nrow(simulation_locations)) {
        punto <- matrix(simulation_locations[point_idx, ], ncol = 2, byrow = T)
        coord_x <- sp::spDistsN1(punto, pt = c(bounds['x', 'min'], punto[[2]]), longlat = !is_projected)
        coord_y <- sp::spDistsN1(punto, pt = c(punto[[1]], bounds['y', 'min']), longlat = !is_projected)
        simulation_grid[point_idx, ] <- c(coord_x, coord_y)
    }

    for(point_idx in 1:nrow(stations_locations)) {
        punto <- matrix(stations_locations[point_idx, ], ncol = 2, byrow = T)
        coord_x <- sp::spDistsN1(punto, pt = c(bounds['x', 'min'], punto[[2]]), longlat = !is_projected)
        coord_y <- sp::spDistsN1(punto, pt = c(punto[[1]], bounds['y', 'min']), longlat = !is_projected)
        stations_grid[point_idx, ] <- c(coord_x, coord_y)
    }

    # simulation_grid <- matrix(simulation_grid, ncol = 2, byrow = T)
    list(simulation_grid = simulation_grid,
         stations_grid = stations_grid)
}


