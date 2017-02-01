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
                                     random_fields_method = 'gaussian',
                                     Rt = NULL,
                                     # grf_method = NULL,
                                     multicombine = F,
                                     always_krig_coefficients = T,
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
    return(list(seasonal_temps_covariates_getter = seasonal_temps_covariates_getter,
                seasonal_prcp_covariates_getter = seasonal_prcp_covariates_getter,
                random_fields_method = rf_function,
                Rt = Rt,
                # grf_method = grf_method,
                always_krig_coefficients = always_krig_coefficients,
                multicombine = multicombine,
                cache_size = cache_size))
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


    #########################################

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

    #########################################

    if(identical(control$random_fields_method, random_field_noise)) {
        cat('Calculating distance matrix for simulation points\n')
        model$simulation_dist_matrix <- as.dist(sp::spDists(simulation_locations), upper = T)
    }

    # Gamma shape and scale
    SH <- c()
    for (station in as.character(model$stations$id)) SH[station] <- model$gamma[[station]]$alpha

    # skip_krigging <- identical(simulation_locations, model$stations) && !control$always_krig_coefficients
    krige_coefficients <- is.null(matching_stations) || control$always_krig_coefficients

    # sim_locations_latlon <- coordinates(sp::spTransform(simulation_locations, CRS("+proj=longlat +datum=WGS84")))
    simulation_locations_grid <- make_distance_grid(simulation_locations)
    simulation_locations <- coordinates(simulation_locations)

    if(!krige_coefficients) {
        coefocc.sim <- matrix(model$coefficients$coefocc[matching_stations, ], nrow = length(matching_stations), byrow = F)
        coefmin.sim <- matrix(model$coefficients$coefmin[matching_stations, ], nrow = length(matching_stations), byrow = F)
        coefmax.sim <- matrix(model$coefficients$coefmax[matching_stations, ], nrow = length(matching_stations), byrow = F)
        SH.sim <- SH[matching_stations]

        colnames(coefocc.sim) <- colnames(model$coefficients$coefocc)
        colnames(coefmin.sim) <- colnames(model$coefficients$coefmin)
        colnames(coefmax.sim) <- colnames(model$coefficients$coefmax)

        rownames(simulation_locations) <- model$stations$id[matching_stations]
        rownames(simulation_locations_grid) <- model$stations$id[matching_stations]
    } else {
        # estimate model coefficients on grid using ordinary kriging (OK) occurrence
        coefocc.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefocc),
                              dimnames = list(NULL, colnames(model$coefficients$coefocc)))
        for (covariate in 1:ncol(coefocc.sim)) {
            coefocc.sim[, covariate] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                                           model$coefficients$coefocc[, covariate]),
                                                                              # sim_locations_latlon))
                                                                              simulation_locations))
        }

        # minimum temperature
        coefmin.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefmin),
                              dimnames = list(NULL, colnames(model$coefficients$coefmin)))
        for (covariate in 1:ncol(coefmin.sim)) {
            coefmin.sim[, covariate] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                                           model$coefficients$coefmin[, covariate]),
                                                                              # sim_locations_latlon))
                                                                              simulation_locations))
        }

        # maximum temperature
        coefmax.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefmax),
                              dimnames = list(NULL, colnames(model$coefficients$coefmax)))
        for (covariate in 1:ncol(coefmax.sim)) {
            coefmax.sim[, covariate] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                                           model$coefficients$coefmax[, covariate]),
                                                                              # sim_locations_latlon))
                                                                              simulation_locations))
        }

        # SH.sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SH), sim_locations_latlon))
        SH.sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SH), simulation_locations))
        # SH.sim <- fields::predict.Krig(fields::Krig(coordinates(coordinates(model$stations)), SH), grid_test)
        if(verbose) cat('Kriged coefficients.\n')

        rownames(simulation_locations) <-sprintf("(%f, %f)", simulation_locations[, 1], simulation_locations[, 2])
        rownames(simulation_locations_grid) <-sprintf("(%f, %f)", simulation_locations[, 1], simulation_locations[, 2])
    }

    # if(identical(control$random_fields_method, cholesky_random_field)) {
    #     cat('Calculate cSigma matrixes.\n')
    #     # cSigma_list <- list()
    #     # # grid     <- list(x = unique(coordinates(simulation_locations_grid)[,1]), y = unique(coordinates(simulation_locations_grid)[,2]))
    #     # # xg       <- fields::make.surface.grid(grid)
    #     #
    #     # model$cSigma <- foreach(month_number=1:12, .multicombine = T) %do% {
    #     #     month_params <- model$month_params[[month_number]]
    #     #     cSigma_list <- list()
    #     #     for(var_name in names(month_params)) {
    #     #
    #     #         # bigSigma <- fields::stationary.cov(xg, theta=month_params[[var_name]][3])
    #     #         bigSigma <- fields::stationary.cov(simulation_locations_grid, theta=month_params[[var_name]][3])
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
    # foreach::getDoParWorkers()

    combine_function <- function(...) {
        invisible(gc())
        abind::abind(..., along = 1)
    }

    simulation_start <- as.Date(start_date) - 1

    # initialize with climatology of the previous day
    # w2 <- suppressWarnings(geoR::grf(nrow(simulation_locations_grid),
    #                                  grid = simulation_locations_grid,
    #                                  cov.model = "exponential",
    #                                  cov.pars = model$month_params[[lubridate::month(simulation_start)]]$prcp[c(2,3)],
    #                                  nugget = model$month_params[[lubridate::month(simulation_start)]]$prcp[1],
    #                                  mean = rep(0, nrow(simulation_locations_grid)), messages = FALSE))
    previous_occ <- (drop(control$random_fields_method(model, simulation_locations_grid, lubridate::month(simulation_start), 'prcp')) > 0) + 0
    # previous_occ <- (noise_generator$generate_noise(lubridate::month(simulation_start), 'prcp') > 0) + 0

    start_climatology <- as.data.frame(model$start_climatology %>% filter(month == lubridate::month(simulation_start), day == lubridate::day(simulation_start)))

    if(!krige_coefficients) {
        previous_tx <- start_climatology[matching_stations, 'tx']
        previous_tn <- start_climatology[matching_stations, 'tn']
    } else {
        previous_tn <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tn), simulation_locations)))
        previous_tx <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tx), simulation_locations)))
    }

    rm(start_climatology)

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

        # temp_covariates <- data.frame(year = years, temp_covariates)
        # prcp_covariates <- data.frame(year = years,
        #                               ST = )

        simulation_dates[season_indexes, paste0("ST", season_number)] <- st_covariates[match(simulation_dates[season_indexes, 'year'], years)]
        simulation_dates[season_indexes, paste0("SMX", season_number)] <- smx_covariates[match(simulation_dates[season_indexes, 'year'], years)]
        simulation_dates[season_indexes, paste0("SMN", season_number)] <- smn_covariates[match(simulation_dates[season_indexes, 'year'], years)]
    }

    # daily_covariates <- as.matrix(t(simulation_dates[, c(5:8, 11:ncol(simulation_dates))]))
    daily_covariates <- as.matrix(t(simulation_dates[, !colnames(simulation_dates) %in% c('date', 'year', 'year_fraction', 'month', 'season')]))

    daily_covariates <- rbind(1, daily_covariates)
    rownames(daily_covariates)[1] <- "(Intercept)"

    SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(coordinates(model$stations))))
    colnames(SC) <- model$stations$id

    for (station in as.character(model$stations$id)) {
        gamma_coef <- model$gamma[[station]]$coef
        SC[, station] <- exp(apply(daily_covariates[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T))/SH[station]
    }

    rm(gamma_coef)

    # simulated_occurrence <- simulated_tx <- simulated_tn <- simulated_prcp <- array(dim = c(nrow(simulation_dates), nrow(simulation_locations_grid)))

    # simulated_climate[1, , , ] <- i
    # simulated_climate[, 1:nrow(simulation_dates), , ] <- 1:nrow(simulation_dates)
    # simulated_climate[, , 1:nrow(simulation_locations_grid), ] <- 1:nrow(simulation_locations_grid)

    # dimnames = list(c('date', 'location')))

    if(verbose) cat('Simulated start climatology.\n')

    simulation_dates <- simulation_dates[, c(1, 2, 3)]

    invisible(gc())



    # noise_generator <- make_noise_proto(model = model,
    #                                     simulation_locations = simulation_locations_grid,
    #                                     control = control,
    #                                     noise_function = control$random_fields_method)

    # gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %dopar% {
    gen_climate <- foreach(i = 1:nsim, .combine = combine_function, .multicombine = control$multicombine, .packages = c('sp')) %dopar% {
    # gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %do% {
        set.seed(realizations_seeds[i])

        noise_generator <- noise_method$new(model = model,
                                            simulation_locations = simulation_locations_grid,
                                            control = control,
                                            noise_function = control$random_fields_method)

        simulated_climate <- array(data = 0.0, dim = c(1, nrow(simulation_dates), nrow(simulation_locations_grid), 3))


        dimnames(simulated_climate)[2] <- list('dates' = format(simulation_dates$date, '%Y-%m-%d'))
        dimnames(simulated_climate)[3] <- list('coordinates' = rownames(simulation_locations))
        dimnames(simulated_climate)[4] <- list('variables' = c('tx', 'tn', 'prcp'))


        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {

            # month_number <- lubridate::month(simulation_dates[d, "date"])
            month_number <- simulation_dates$month[d]


            simulation_matrix <- cbind(previous_occ, previous_tn, previous_tx, matrix(daily_covariates[, d], ncol = nrow(daily_covariates), nrow = length(previous_tx), byrow = T))
            colnames(simulation_matrix) <- c("prcp_occ_prev", "tn_prev", "tx_prev", rownames(daily_covariates))


            # mu.occ <- apply(simulation_matrix[, colnames(coefocc.sim)] * coefocc.sim, 1, sum, na.rm = T)
            mu_occ <- rowSums(simulation_matrix[, colnames(coefocc.sim)] * coefocc.sim)
            # occ_noise <- control$random_fields_method(model, simulation_locations_grid, month_number, 'prcp')
            occ_noise <- noise_generator$generate_noise(month_number, 'prcp')
            current_occ <- ((mu_occ + occ_noise) > 0) + 0

            # simulated_climate[1, d, , 'prcp'] <- ((mu_occ + occ_noise) > 0) + 0

            simulation_matrix <- cbind(prcp_occ = current_occ, simulation_matrix)

            # mu.min <- apply(simulation_matrix[, colnames(coefmin.sim)] * coefmin.sim, 1, sum, na.rm = T)
            # mu.max <- apply(simulation_matrix[, colnames(coefmax.sim)] * coefmax.sim, 1, sum, na.rm = T)
            mu_tn <- rowSums(simulation_matrix[, colnames(coefmin.sim)] * coefmin.sim)
            # tn_noise <- control$random_fields_method(model, simulation_locations_grid, month_number, 'tn')
            tn_noise <- noise_generator$generate_noise(month_number, 'tn')
            current_tn <- signif(mu_tn + tn_noise, digits = 4)

            mu_tx <- rowSums(simulation_matrix[, colnames(coefmax.sim)] * coefmax.sim)
            # tx_noise <- control$random_fields_method(model, simulation_locations_grid, month_number, 'tx')
            tx_noise <- noise_generator$generate_noise(month_number, 'tx')
            current_tx <- signif(mu_tx + tx_noise, digits = 4)

            # simulated_tx[d, ] <- signif(mu.max + w.max, digits = 4)

            daily_retries <- 0
            # 1/10000
            # while (min(simulated_climate[1, d, , 'tx'] - simulated_climate[1, d, , 'tn']) < 0.5) {
            while (min(current_tx - current_tn) < 0.5 && daily_retries < 100) {
                temps_retries <- temps_retries + 1

                # tn_noise <- control$random_fields_method(model, simulation_locations_grid, month_number, 'tn')
                # tx_noise <- control$random_fields_method(model, simulation_locations_grid, month_number, 'tx')
                tn_noise <- noise_generator$generate_noise(month_number, 'tn')
                tx_noise <- noise_generator$generate_noise(month_number, 'tx')

                current_tx <- signif(mu_tx + tx_noise, digits = 4)
                current_tn <- signif(mu_tn + tn_noise, digits = 4)

                daily_retries <- daily_retries + 1
                # simulated_climate[1, d, , 'tn'] <- signif(mu_tn + tn_noise, digits = 4)
                # simulated_climate[1, d, , 'tx'] <- signif(mu_tx + tx_noise, digits = 4)
            }
            if(daily_retries >= 100) stop('Failed to simulate random noise that doesn\'t violate the constraint of max. t. > min. t')

            simulated_climate[1, d, , 'prcp'] <- current_occ
            simulated_climate[1, d, , 'tn'] <- current_tn
            simulated_climate[1, d, , 'tx'] <- current_tx

            # SIMocc.old <- simulated_occurrence[d, ]
            # SIMmax.old <- simulated_tx[d, ]
            # SIMmin.old <- simulated_tn[d, ]
            previous_occ <- current_occ
            previous_tx <- current_tx
            previous_tn <- current_tn

            ## amounts
            # SC.sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SC[d, ]), simulation_locations))
            #
            # # w3 <- suppressWarnings(geoR::grf(nrow(simulation_locations_grid), grid = simulation_locations_grid, cov.model = "exponential", cov.pars = c(p$prcp[2], p$prcp[3]), nugget = p$prcp[1], mean = rep(0,
            # #     nrow(simulation_locations_grid)), messages = FALSE)$data)
            #
            #
            # rainy_locations <- simulated_climate[1, d, , 'prcp'] > 0
            #
            # if(any(rainy_locations)) {
            #     w3 <- control$random_fields_method(model, simulation_locations_grid, month_number, 'prcp')
            #     simulated_climate[1, d, rainy_locations, 'prcp'] <- signif(qgamma(pnorm(w3), shape = SH.sim, scale = SC.sim), digits = 4)[rainy_locations]
            # }
            # simulated_prcp[d, ] <- signif(qgamma(pnorm(w3), shape = SH.sim, scale = SC.sim), digits = 4)
            if(verbose) cat(paste0("\r", i, ": ", d, "/", ncol(daily_covariates), ". Retries: ", temps_retries, '       '))

            if(d %% 1000 == 0) invisible(gc())
        }

        simulated_climate[1, , , 'prcp'][simulated_climate[1, , , 'prcp'] < model$control$prcp_occurrence_threshold] <- 0

        simulated_climate
    }


    attr(gen_climate, 'krigged_coefficients') <- list(
        'tx' = coefmax.sim,
        'tn' = coefmin.sim,
        'occ' = coefocc.sim
    )

    attr(gen_climate, 'seasonal_covariates') <- list(
        'tx' = smx_covariates,
        'tn' = smn_covariates,
        'prcp' = st_covariates
    )

    attr(gen_climate, 'realizations_seeds') <- realizations_seeds

    attr(gen_climate, 'distance_grid') <- simulation_locations_grid

    # save(gen_climate, file='gen_climate.RData')


    # if (nsim == 1) {
    #     # gen_climate <- list(gen_climate)
    #     return(gen_climate[1, , , ])
    # }

    gen_climate
}
