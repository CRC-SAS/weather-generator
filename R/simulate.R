partial <- function(f, ...) {
    l <- list(...)
    function(...) {
        do.call(f, c(list(...), l))
    }
}

# get_seasonal_covariate <- function(year, season_number, season_values, quantile_probs = c(0, 0.33, 0.66, 1), levels_probabilities = c(1/3, 1/3, 1/3)) {
#     splitted <- cut(season_values, breaks = quantile(season_values, probs = quantile_probs), include.lowest = T)
#     values_level <- season_values[splitted == sample(levels(splitted), 1, prob = levels_probabilities)]
#     if(length(values_level) == 1) return(values_level)
#     return(base::sample(values_level, 1))
# }

get_seasonal_covariate <- function(years, season_number, seasonal_values, variable) {
    seasonal_values <- seasonal_values %>% filter(season == season_number) %>% select(year, value = ends_with(variable))
    if(all(years %in% seasonal_values$year)) {
        return(seasonal_values[match(years, seasonal_values$year), 'value'])
    }
    warning('Seasonal values not found for some years, starting series at the first year')
    return(seasonal_values$value[1:length(years)])
}

get_rainfall_seasonal_covariate <- function(years, season_number, seasonal_values) {
    get_seasonal_covariate(years, season_number, seasonal_values, 'prcp')
}

# get_temperatures_seasonal_covariate <- function(year, season_number, season_tx_values, season_tn_values, quantile_probs = c(0, 0.33, 0.66, 1),
#                                                levels_probabilities = c(0.25, 0.35, 0.4)) {
#     splitted_tx <- cut(season_tx_values, breaks = quantile(season_tx_values, probs = quantile_probs), include.lowest = T)
#     splitted_tn <- cut(season_tn_values, breaks = quantile(season_tn_values, probs = quantile_probs), include.lowest = T)
#
#     values_level_index <- sample(seq_along(levels(splitted_tx)), 1, prob = levels_probabilities)
#     tx_sampled_level <- season_tx_values[splitted_tx == levels(splitted_tx)[values_level_index]]
#     tn_sampled_level <- season_tn_values[splitted_tn == levels(splitted_tn)[values_level_index]]
#     return(list(tx = ifelse(length(tx_sampled_level) > 1, sample(tx_sampled_level, 1), tx_sampled_level),
#                 tn = ifelse(length(tn_sampled_level) > 1, sample(tn_sampled_level, 1), tn_sampled_level)))
# }
get_temperatures_seasonal_covariate <- function(years, season_number, seasonal_values) {
    return(list(tx = get_seasonal_covariate(years, season_number, seasonal_values, 'tx'),
                tn = get_seasonal_covariate(years, season_number, seasonal_values, 'tn')))
}

#' @title Simulations control configuration
#' @description Provides fine control of different parameters that will be used to create new weather series.
#' @export
glmwgen_simulation_control <- function(seasonal_temps_covariates_getter = get_temperatures_seasonal_covariate,
                                       seasonal_prcp_covariates_getter = get_rainfall_seasonal_covariate,
                                       random_fields_method = 'rf',
                                       minimum_temperatures_difference_threshold = 0.1,
                                       Rt = NULL,
                                       # grf_method = NULL,
                                       multicombine = F,
                                       always_krig_coefficients = T,
                                       interpolation_method = c('idw', 'kriging'),
                                       cache_size = 1) {
    rf_function <- NULL
    if(!is.function(random_fields_method)) {
        # if(startsWith(random_fields_method, 'gauss')) rf_function <- gaussian_random_field
        if(startsWith(random_fields_method, 'rf')) rf_function <- random_field_noise
        if(startsWith(random_fields_method, 'chol')) rf_function <- cholesky_random_field
        if(startsWith(random_fields_method, 'rn')) rf_function <- rnorm_noise
        if(startsWith(random_fields_method, 'mv')) rf_function <- mvrnorm_noise
    } else {
        rf_function <- random_fields_method
    }

    interpolation_method <- match.arg(interpolation_method)

    interpolation_method <- c('kriging' = krige_covariate_automap, 'idw' = idw_covariate)[[interpolation_method]]

    return(list(seasonal_temps_covariates_getter = seasonal_temps_covariates_getter,
                seasonal_prcp_covariates_getter = seasonal_prcp_covariates_getter,
                random_fields_method = rf_function,
                Rt = Rt,
                # grf_method = grf_method,
                always_krig_coefficients = always_krig_coefficients,
                multicombine = multicombine,
                interpolation_method = interpolation_method,
                cache_size = cache_size,
                minimum_temperatures_difference_threshold = minimum_temperatures_difference_threshold))
}

#' @title Simulates new weather trajectories
#' @description Simulates new weather trajectories.
#' @param object A glmwgen model.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#'
#'          Either NULL or an integer that will be used in a call to set.seed before simulating the response vectors.
#'          If set, the value is saved as the "seed" attribute of the returned value. The default, NULL will not change the random generator state.
#' @param start_date a start date in text format (will be converted using as.Date) or a date object.
#' @param end_date an end date in text format (will be converted using as.Date) or a date object.
#' @param simulation_locations a SpatialPoints object with the points at which weather should be simulated.
#'          If not set, the locations used to fit the model will be used.
#' @param control a glmwgen simulation control list.
#' @import dplyr
#' @import foreach
#' @export
simulate.glmwgen <- function(object, nsim = 1, seed = NULL, start_date = NA, end_date = NA, simulation_locations = NULL, control = glmwgen_simulation_control(),
                             verbose = T) {
    model <- object

    if(class(object) != 'glmwgen') {
        stop(paste('Received a model of class', class(object), 'and a model of class "glmwgen" was expected.'))
    }

    if(!is.null(seed)) set.seed(seed)

    if (is.null(simulation_locations)) {
        simulation_locations <- model$stations

        if(verbose) cat('Using the coordinates of the fitted stations as simulation locations\n')
    }

    realizations_seeds <- ceiling(runif(min = 1, max = 10000000, n = nsim))


    if(!('proj4string' %in% names(attributes(simulation_locations)))) {
        warning('simulation_locations is not a spacial points object, attempting to convert it with stations projection string.')
        simulation_locations <- SpatialPoints(simulation_locations, proj4string=model$stations_proj4string)
    }

    # Check if all the simulation locations match a station.
    distance_to_stations <- round(sp::spDists(simulation_locations, model$stations), 1)
    matching_stations <- NULL

    if(all(apply(distance_to_stations, 1, min) < 0.5)) {
        matching_stations <- apply(distance_to_stations, 1, which.min)
    }

    # Save info about wether the CRS is a projected one or not.
    is_projected <- is.projected(simulation_locations)

    #########################################

    if(end_date <= start_date) {
        stop('End date should be greater than start date')
    }

    if(nsim < 1) stop('Number of simulations should be greater than one')

    simulation_dates <- data.frame(date = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "days"))

    if (is.null(control$Rt)) {
        # TODO: shouldn't this be = 1 once we start simulating?
        # Rt <- seq(from = -1, to = 1, length.out = nrow(simulation_dates))
        Rt <- 0
    } else {
        Rt <- control$Rt
        if (length(Rt) == 1)
            Rt <- rep(Rt, times = nrow(simulation_dates))
        if (length(Rt) != nrow(simulation_dates))
            stop("The specified Rt parameter length differs from the simulation dates series length.")
    }

    simulation_dates <- simulation_dates %>%
        mutate(year = lubridate::year(date),
               month = lubridate::month(date),
               year_fraction = 2 * pi * lubridate::yday(date) / ifelse(lubridate::leap_year(date), 366, 365),
               ct = cos(year_fraction),
               st = sin(year_fraction),
               Rt = Rt,
               season = ceiling(lubridate::month(date)/3),
               ST1 = 0, ST2 = 0, ST3 = 0, ST4 = 0,
               SMX1 = 0, SMX2 = 0, SMX3 = 0, SMX4 = 0,
               SMN1 = 0, SMN2 = 0, SMN3 = 0, SMN4 = 0)

    rm(Rt)

    for (season_number in unique(simulation_dates$season)) {
        season_indexes <- simulation_dates$season == season_number
        # season_values <- d_seatot[d_seatot$season == season_number, ]
        season_tx <- unique(model$seasonal$tx[[season_number]])
        season_tn <- unique(model$seasonal$tn[[season_number]])
        season_prcp <- unique(model$seasonal$prcp[[season_number]])

        years <- unique(simulation_dates$year[season_indexes])

        temp_covariates <- control$seasonal_temps_covariates_getter(years, season_number, model$seasonal)
        smx_covariates <- temp_covariates[['tx']]
        smn_covariates <- temp_covariates[['tn']]
        st_covariates <- control$seasonal_prcp_covariates_getter(years, season_number, model$seasonal)

        simulation_dates[season_indexes, paste0("ST", season_number)] <- st_covariates[match(simulation_dates[season_indexes, 'year'], years)]
        simulation_dates[season_indexes, paste0("SMX", season_number)] <- smx_covariates[match(simulation_dates[season_indexes, 'year'], years)]
        simulation_dates[season_indexes, paste0("SMN", season_number)] <- smn_covariates[match(simulation_dates[season_indexes, 'year'], years)]
    }

    daily_covariates <- as.matrix(t(simulation_dates[, !colnames(simulation_dates) %in% c('date', 'year', 'year_fraction', 'month', 'season')]))

    daily_covariates <- rbind(1, daily_covariates)
    rownames(daily_covariates)[1] <- "(Intercept)"

    #########################################

    if(identical(control$random_fields_method, glmwgen:::random_field_noise)) {
        cat('Calculating distance matrix for simulation points\n')
        model$simulation_dist_matrix <- as.dist(sp::spDists(simulation_locations), upper = T)
    }

    if(identical(control$interpolation_method, glmwgen:::idw_covariate)) {
        model$idw_weights <- 1 / (sp::spDists(simulation_locations, model$stations)) ^ 2
        model$idw_weights[is.infinite(model$idw_weights)] <- 1
        model$idw_weights <- model$idw_weights / rowSums(model$idw_weights)
    }

    # Gamma shape
    SH <- c()
    for (station in as.character(model$stations$id)) SH[station] <- model$gamma[[station]]$alpha

    # Check wether we should krige coefficients or not.
    krige_coefficients <- is.null(matching_stations) || control$always_krig_coefficients

    # Save original coordinates.
    simulation_coordinates <- coordinates(simulation_locations)

    # Create a gridded representation of the simulation points and fitted stations.
    # This is similar to a Gauss-Kruger projection.
    projections_grid <- make_distance_grid(model$stations, simulation_locations)
    simulation_coordinates_grid <- projections_grid$simulation_grid

    # stations_in_grid <- apply(spDists(model$stations, simulation_locations), 1, which.min)
    # stations_grid_distance <- apply(spDists(model$stations, simulation_locations), 1, min)
    # station_coordinates_grid <- simulation_coordinates_grid[stations_in_grid, ]
    station_coordinates_grid <- projections_grid$stations_grid

    if(!krige_coefficients) {
        coefocc_sim <- matrix(model$coefficients$coefocc[matching_stations, ], nrow = length(matching_stations), byrow = F)
        coefmin_sim <- matrix(model$coefficients$coefmin[matching_stations, ], nrow = length(matching_stations), byrow = F)
        coefmax_sim <- matrix(model$coefficients$coefmax[matching_stations, ], nrow = length(matching_stations), byrow = F)
        SH_sim <- SH[matching_stations]

        colnames(coefocc_sim) <- colnames(model$coefficients$coefocc)
        colnames(coefmin_sim) <- colnames(model$coefficients$coefmin)
        colnames(coefmax_sim) <- colnames(model$coefficients$coefmax)

        rownames(simulation_coordinates) <- model$stations$id[matching_stations]
        rownames(simulation_coordinates_grid) <- model$stations$id[matching_stations]
    } else {
        stations_krige_sp <- model$stations
        simulation_krige_sp <- sp::SpatialPointsDataFrame(coordinates(simulation_locations),
                                                          data = data.frame(loc_id = 1:nrow(simulation_locations)),
                                                          proj4string = CRS(proj4string(simulation_locations)))

        if(!is.projected(stations_krige_sp)) {
            projection_string <- "+proj=tpeqd +lat_1=%f +lon_1=%f +lat_2=%f +lon_2=%f +x_0=%f +y_0=%f +ellps=intl +units=m +no_defs"
            # projection_string <- "+proj=tmerc +lat_0=%f +lon_0=%f +x_0=%f +y_0=%f +ellps=intl +datum=WGS84 +units=m +no_defs"
            # Get locations bounds.
            bounds <- sp::bbox(rbind(coordinates(stations_krige_sp),
                                     coordinates(simulation_krige_sp)))
            # Format projection string and make it a CRS.
            projection_crs <- sp::CRS(sprintf(projection_string, bounds[2, 'min'], bounds[1, 'min'], bounds[2, 'max'], bounds[1, 'max'], 0, 0))
            # projection_crs <- sp::CRS(sprintf(projection_string, bounds[2, 'min'], bounds[1, 'min'], 0, 0))
            .stations_krige_sp <- sp::spTransform(stations_krige_sp, projection_crs)
            .simulation_krige_sp <- sp::spTransform(simulation_krige_sp, projection_crs)

            coords_offset <- max(abs(apply(rbind(coordinates(.simulation_krige_sp), coordinates(.stations_krige_sp)), 2, min))) + 1000

            projection_crs <- sp::CRS(sprintf(projection_string, bounds[2, 'min'], bounds[1, 'min'], bounds[2, 'max'], bounds[1, 'max'], ceiling(coords_offset), ceiling(coords_offset)))
            # projection_crs <- sp::CRS(sprintf(projection_string, bounds[2, 'min'], bounds[1, 'min'], min_coords[1], min_coords[2]))
            stations_krige_sp <- sp::spTransform(stations_krige_sp, projection_crs)
            simulation_krige_sp <- sp::spTransform(simulation_krige_sp, projection_crs)

            rm(.stations_krige_sp, .simulation_krige_sp)
            cat('Projected stations and simulation locations to interpolate coefficients\n')
        }

        coefmin_sim <- apply(model$coefficients$coefmin, 2, control$interpolation_method, model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp)
        coefocc_sim <- apply(model$coefficients$coefocc, 2, control$interpolation_method, model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp)
        coefmax_sim <- apply(model$coefficients$coefmax, 2, control$interpolation_method, model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp)
        cat('Interpolating gamma shape, this may take a while...')

        SH_sim <- control$interpolation_method(model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp, SH)

        SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(coordinates(model$stations))))
        colnames(SC) <- model$stations$id

        for (station in as.character(model$stations$id)) {
            gamma_coef <- model$gamma[[station]]$coef
            SC[, station] <- exp(apply(daily_covariates[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T)) / SH[station]
        }

        # Andrew, new
        start_time <- Sys.time()

        `%op%` <- if(identical(control$interpolation_method, glmwgen:::idw_covariate)) `%do%` else `%dopar%`

        SC_sim <- foreach(day_idx = 1:nrow(SC), .combine = rbind, .packages = c('sp')) %op% {
            control$interpolation_method(model = model,
                                         stations_locations = stations_krige_sp,
                                         simulation_locations = simulation_krige_sp,
                                         SC[day_idx, ])
        }
        cat(sprintf('done (took %.2f minutes)\n', as.numeric(difftime(Sys.time(), start_time, units = 'min'))))

        # Fede
        # coefsc <- sapply(as.character(stations_krige_sp$id),
        #                  function(station_id) model$gamma[[station_id]]$coef)
        # coefsc_sim <- apply(coefsc, 1, control$interpolation_method,
        #                     model = model, stations_locations = stations_krige_sp,
        #                     simulation_locations = simulation_krige_sp)
        #
        # SC_sim <- t(exp(coefsc_sim %*% daily_covariates[names(gamma_coef), ]) / SH_sim)

        # Andrew, old
        # SC_sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SC[d, ]), simulation_coordinates))

        rm(gamma_coef)

        if(verbose) cat('Interpolated coefficients.\n')

        if(!is.null(matching_stations)) {
            rownames(simulation_coordinates) <- model$stations$id[matching_stations]
            rownames(simulation_coordinates_grid) <- model$stations$id[matching_stations]
        } else {
            rownames(simulation_coordinates) <- sprintf("(%f, %f)", simulation_coordinates[, 1], simulation_coordinates[, 2])
            rownames(simulation_coordinates_grid) <- sprintf("(%f, %f)", simulation_coordinates[, 1], simulation_coordinates[, 2])
        }
        rownames(coefmin_sim) <- rownames(coefmax_sim) <- rownames(coefocc_sim) <- colnames(SC_sim) <- rownames(simulation_coordinates_grid)
    }

    # if(identical(control$random_fields_method, cholesky_random_field)) {
    #     cat('Calculate cSigma matrixes.\n')
    #     # cSigma_list <- list()
    #     # # grid     <- list(x = unique(coordinates(simulation_coordinates_grid)[,1]), y = unique(coordinates(simulation_coordinates_grid)[,2]))
    #     # # xg       <- fields::make.surface.grid(grid)
    #     #
    #     # model$cSigma <- foreach(month_number=1:12, .multicombine = T) %do% {
    #     #     month_params <- model$month_params[[month_number]]
    #     #     cSigma_list <- list()
    #     #     for(var_name in names(month_params)) {
    #     #
    #     #         # bigSigma <- fields::stationary.cov(xg, theta=month_params[[var_name]][3])
    #     #         bigSigma <- fields::stationary.cov(simulation_coordinates_grid, theta=month_params[[var_name]][3])
    #     #         cSigma_list[[var_name]] <- t(chol(bigSigma))
    #     #     }
    #     #     print(paste0('Chol transform. Month = ', month, '. Variable = ', var_name))
    #     #
    #     #     cSigma_list
    #     # }
    #     #
    #     # # names(model$cSigma) <- unique(lubridate::month(simulation_dates$date))
    #     #
    #     # # rm(grid, xg, bigSigma)
    #     # rm(bigSigma)
    # }

    # Register a sequential backend if the user didn't register a parallel
    # in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        foreach::registerDoSEQ()
    }
    process_pid <- Sys.getpid()

    combine_function <- function(...) {
        invisible(gc())
        abind::abind(..., along = 1)
    }

    simulation_start <- as.Date(start_date) - 1

    previous_occ <- (drop(control$random_fields_method(model, simulation_coordinates_grid, lubridate::month(simulation_start), 'prcp')) > 0) + 0

    start_climatology <- as.data.frame(model$start_climatology %>% filter(month == lubridate::month(simulation_start), day == lubridate::day(simulation_start)))

    if(!krige_coefficients) {
        previous_tx <- start_climatology[matching_stations, 'tx']
        previous_tn <- start_climatology[matching_stations, 'tn']
    } else {
        # previous_tn <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tn), simulation_coordinates)))
        # previous_tx <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tx), simulation_coordinates)))
        previous_tn <- control$interpolation_method(model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp, start_climatology$tn)
        previous_tx <- control$interpolation_method(model = model, stations_locations = stations_krige_sp, simulation_locations = simulation_krige_sp, start_climatology$tx)
    }

    rm(start_climatology)



    if(verbose) cat('Simulated start climatology.\n')

    simulation_dates <- simulation_dates[, c(1, 2, 3)]

    invisible(gc())

    # gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %dopar% {
    gen_climate <- foreach(i = 1:nsim, .combine = combine_function, .multicombine = control$multicombine, .packages = c('sp')) %dopar% {
        set.seed(realizations_seeds[i])

        noise_generator <- noise_method$new(model = model,
                                            simulation_locations = simulation_coordinates_grid,
                                            control = control,
                                            noise_function = control$random_fields_method)

        simulated_climate <- array(data = 0.0, dim = c(1, nrow(simulation_dates), nrow(simulation_coordinates_grid), 3))


        dimnames(simulated_climate)[2] <- list('dates' = format(simulation_dates$date, '%Y-%m-%d'))
        dimnames(simulated_climate)[3] <- list('coordinates' = rownames(simulation_coordinates))
        dimnames(simulated_climate)[4] <- list('variables' = c('tx', 'tn', 'prcp'))


        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {
            month_number <- simulation_dates$month[d]


            simulation_matrix <- cbind(previous_occ, previous_tn, previous_tx, matrix(daily_covariates[, d], ncol = nrow(daily_covariates), nrow = length(previous_tx), byrow = T))
            colnames(simulation_matrix) <- c("prcp_occ_prev", "tn_prev", "tx_prev", rownames(daily_covariates))


            mu_occ <- rowSums(simulation_matrix[, colnames(coefocc_sim)] * coefocc_sim)
            occ_noise <- noise_generator$generate_noise(month_number, 'prcp')
            simulated_climate[1, d, , 'prcp'] <- ((mu_occ + occ_noise) > 0) + 0


            simulation_matrix <- cbind(prcp_occ = simulated_climate[1, d, , 'prcp'], simulation_matrix)


            mu_tn <- rowSums(simulation_matrix[, colnames(coefmin_sim)] * coefmin_sim)
            tn_noise <- noise_generator$generate_noise(month_number, 'tn')
            simulated_climate[1, d, , 'tn'] <- signif(mu_tn + tn_noise, digits = 4)

            mu_tx <- rowSums(simulation_matrix[, colnames(coefmax_sim)] * coefmax_sim)
            tx_noise <- noise_generator$generate_noise(month_number, 'tx')
            simulated_climate[1, d, , 'tx'] <- signif(mu_tx + tx_noise, digits = 4)


            daily_retries <- 0
            # retries_array <- array(NA, dim = c(100, length(simulated_climate[1, d, , 'tx'])))
            while (min(simulated_climate[1, d, , 'tx'] - simulated_climate[1, d, , 'tn']) < control$minimum_temperatures_difference_threshold && daily_retries < 100) {
                temps_retries <- temps_retries + 1
                daily_retries <- daily_retries + 1

                # retries_array[daily_retries, ] <- simulated_climate[1, d, , 'tx'] - simulated_climate[1, d, , 'tn']

                tn_noise <- noise_generator$generate_noise(month_number, 'tn')
                tx_noise <- noise_generator$generate_noise(month_number, 'tx')

                simulated_climate[1, d, , 'tx'] <- signif(mu_tx + tx_noise, digits = 4)
                simulated_climate[1, d, , 'tn'] <- signif(mu_tn + tn_noise, digits = 4)
            }
            if(daily_retries >= 100) {
                cat('Failed to simulate random noise that doesn\'t violate the constraint of max. temp. > min. temp.')
                # temps_diff <- apply(retries_array, 2, function(x) x < control$minimum_temperatures_difference_threshold)
                # temps_diff <- apply(temps_diff, 2, sum)
                # fields::quilt.plot(simulation_coordinates, temps_diff)
                # text(coordinates(model$stations), label = model$stations$id, col = 'white')
                # title(main = sprintf('Day %03d', d))
                return(NULL)
            }

            previous_occ <- simulated_climate[1, d, , 'prcp']
            previous_tx <- simulated_climate[1, d, , 'tx']
            previous_tn <- simulated_climate[1, d, , 'tn']


            rainy_locations <- simulated_climate[1, d, , 'prcp'] > 0

            if(any(rainy_locations)) {
                # amounts

                # w3 <- suppressWarnings(geoR::grf(nrow(simulation_coordinates_grid), grid = simulation_coordinates_grid, cov.model = "exponential", cov.pars = c(p$prcp[2], p$prcp[3]), nugget = p$prcp[1], mean = rep(0,
                #     nrow(simulation_coordinates_grid)), messages = FALSE)$data)
                #
                # w3 <- control$random_fields_method(model, simulation_coordinates_grid, month_number, 'prcp')
                amt_noise <- noise_generator$generate_noise(month_number, 'prcp')

                simulated_climate[1, d, rainy_locations, 'prcp'] <- signif(qgamma(pnorm(amt_noise), shape = SH_sim, scale = SC_sim[d, ]), digits = 4)[rainy_locations]
                # simulated_prcp[d, ] <- signif(qgamma(pnorm(w3), shape = SH_sim, scale = SC.sim), digits = 4)
            }
            if(verbose && d %% 10 == 0) cat(paste0("\r", i, ": ", d, "/", ncol(daily_covariates), ". Retries: ", temps_retries, '       '))

            if(d %% 30 == 0) invisible(gc())
        }

        simulated_climate[1, , , 'prcp'][simulated_climate[1, , , 'prcp'] < model$control$prcp_occurrence_threshold] <- 0

        simulated_climate
    }


    attr(gen_climate, 'krigged_coefficients') <- list(
        'tx' = coefmax_sim,
        'tn' = coefmin_sim,
        'occ' = coefocc_sim
    )

    attr(gen_climate, 'seasonal_covariates') <- list(
        'tx' = smx_covariates,
        'tn' = smn_covariates,
        'prcp' = st_covariates
    )

    attr(gen_climate, 'realizations_seeds') <- realizations_seeds
    attr(gen_climate, 'simulation_coordinates') <- simulation_coordinates

    if('simulation_dist_matrix' %in% names(model)) {
        model[['simulation_dist_matrix']] <- NULL
    }

    attr(gen_climate, 'model') <- model

    class(gen_climate) <- c(class(gen_climate), 'glmwgen.climate')

    gen_climate
}
