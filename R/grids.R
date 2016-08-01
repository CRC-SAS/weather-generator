make_distance_grid <- function(simulation_locations) {
    is_projected <- is.projected(simulation_locations)
    simulation_locations <- coordinates(simulation_locations)
    bounds <- sp::bbox(simulation_locations)
    dist_grid <- matrix(nrow = nrow(simulation_locations), ncol = 2)

    for(point_idx in 1:nrow(simulation_locations)) {
        punto <- matrix(simulation_locations[point_idx, ], ncol = 2, byrow = T)
        coord_x <- sp::spDistsN1(punto, pt = c(bounds['x', 'min'], punto[[2]]), longlat = !is_projected)
        coord_y <- sp::spDistsN1(punto, pt = c(punto[[1]], bounds['y', 'min']), longlat = !is_projected)
        dist_grid[point_idx, ] <- c(coord_x, coord_y)
    }

    # dist_grid <- matrix(dist_grid, ncol = 2, byrow = T)
    dist_grid
}


