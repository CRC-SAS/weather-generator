
#' @title Simulates new weather trajectories in stations
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
sim.stns.glmwgen <- function(object, nsim = 1, seed = NULL, start_date = NA, end_date = NA,
                             control = glmwgen:::glmwgen_simulation_control(),  verbose = T) {
    model <- object

    if(class(object) != 'glmwgen') {
        stop(paste('Received a model of class', class(object), 'and a model of class "glmwgen" was expected.'))
    }

    if(!is.null(seed)) set.seed(seed)

    simulation_locations <- model$stations

    realizations_seeds <- ceiling(runif(min = 1, max = 10000000, n = nsim))


    if(!('proj4string' %in% names(attributes(simulation_locations)))) {
        warning('simulation_locations is not a spacial points object, attempting to convert it with stations projection string.')
        simulation_locations <- SpatialPoints(simulation_locations, proj4string=model$stations_proj4string)
    }

    # Save info about wether the CRS is a projected one or not.
    is_projected <- sp::is.projected(simulation_locations)

    #########################################

    if(end_date <= start_date) {
        stop('End date should be greater than start date')
    }

    if(nsim < 1) stop('Number of simulations should be greater than one')

    if(identical(control$random_fields_method, glmwgen:::random_field_noise)) {
        if(verbose) cat("There isn't possible to use random_field_noise function to calculate previous_occ. Switching to rnorm_noise function!\n")
        control$random_fields_method <- glmwgen:::rnorm_noise
    }

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
        # season_tx <- unique(model$seasonal$tx[[season_number]])
        # season_tn <- unique(model$seasonal$tn[[season_number]])
        # season_prcp <- unique(model$seasonal$prcp[[season_number]])

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

    # Gamma shape
    SH <- c()
    for (station in as.character(model$stations$id)) SH[station] <- model$gamma[[station]]$alpha

    stations <- seq_len(nrow(simulation_locations))  # en lugar de matching_stations

    # Save original coordinates.
    simulation_coordinates <- sp::coordinates(simulation_locations)
    rownames(simulation_coordinates) <- model$stations$id[stations]


        coefocc_sim <- matrix(model$coefficients$coefocc[stations, ], nrow = length(stations), byrow = F)
        coefmin_sim <- matrix(model$coefficients$coefmin[stations, ], nrow = length(stations), byrow = F)
        coefmax_sim <- matrix(model$coefficients$coefmax[stations, ], nrow = length(stations), byrow = F)
        SH <- SH[stations]

        colnames(coefocc_sim) <- colnames(model$coefficients$coefocc)
        colnames(coefmin_sim) <- colnames(model$coefficients$coefmin)
        colnames(coefmax_sim) <- colnames(model$coefficients$coefmax)

        SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(sp::coordinates(model$stations))))
        colnames(SC) <- model$stations$id

        for (station in as.character(model$stations$id)) {
            gamma_coef <- model$gamma[[station]]$coef
            SC[, station] <- exp(apply(daily_covariates[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T)) / SH[station]
        }

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

    previous_occ <- (drop(control$random_fields_method(model, simulation_coordinates, lubridate::month(simulation_start), 'prcp')) > 0) + 0

    start_climatology <- as.data.frame(model$start_climatology %>% filter(month == lubridate::month(simulation_start), day == lubridate::day(simulation_start)))

    stations <- seq_len(nrow(simulation_locations))
    previous_tx <- start_climatology[stations, 'tx']
    previous_tn <- start_climatology[stations, 'tn']

    rm(start_climatology)



    if(verbose) cat('Simulated start climatology.\n')

    simulation_dates <- simulation_dates[, c(1, 2, 3)]

    invisible(gc())

    # gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %dopar% {
    gen_climate <- foreach(i = 1:nsim, .combine = combine_function, .multicombine = control$multicombine, .packages = c('sp')) %dopar% {
        set.seed(realizations_seeds[i])

        noise_generator <- glmwgen:::noise_method$new(model = model,
                                                      simulation_locations = simulation_coordinates,
                                                      control = control,
                                                      noise_function = control$random_fields_method)

        simulated_climate <- array(data = 0.0, dim = c(1, nrow(simulation_dates), nrow(simulation_coordinates), 3))


        dimnames(simulated_climate)[2] <- list('dates' = format(simulation_dates$date, '%Y-%m-%d'))
        dimnames(simulated_climate)[3] <- list('coordinates' = rownames(simulation_coordinates))
        dimnames(simulated_climate)[4] <- list('variables' = c('tx', 'tn', 'prcp'))


        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {
            month_number <- simulation_dates$month[d]


            simulation_matrix <- cbind(previous_occ, previous_tn, previous_tx,
                                       matrix(daily_covariates[, d], ncol = nrow(daily_covariates), nrow = length(previous_tx), byrow = T))
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
                if(verbose) cat('Failed to simulate random noise that doesn\'t violate the constraint of max. temp. > min. temp.')
                # temps_diff <- apply(retries_array, 2, function(x) x < control$minimum_temperatures_difference_threshold)
                # temps_diff <- apply(temps_diff, 2, sum)
                # fields::quilt.plot(simulation_coordinates, temps_diff)
                # text(sp::coordinates(model$stations), label = model$stations$id, col = 'white')
                # title(main = sprintf('Day %03d', d))
                return(NULL)
            }

            # aca considerar los lags!!
            previous_occ <- simulated_climate[1, d, , 'prcp']
            previous_tx <- simulated_climate[1, d, , 'tx']
            previous_tn <- simulated_climate[1, d, , 'tn']


            rainy_locations <- simulated_climate[1, d, , 'prcp'] > 0

            if(any(rainy_locations)) {
                # amounts generation

                # w3 <- suppressWarnings(geoR::grf(nrow(simulation_coordinates_grid), grid = simulation_coordinates_grid, cov.model = "exponential", cov.pars = c(p$prcp[2], p$prcp[3]), nugget = p$prcp[1], mean = rep(0,
                #     nrow(simulation_coordinates_grid)), messages = FALSE)$data)
                #
                # w3 <- control$random_fields_method(model, simulation_coordinates_grid, month_number, 'prcp')
                amt_noise <- noise_generator$generate_noise(month_number, 'prcp')

                simulated_climate[1, d, rainy_locations, 'prcp'] <- signif(qgamma(pnorm(amt_noise), shape = SH, scale = SC[d, ]), digits = 4)[rainy_locations]
                # simulated_prcp[d, ] <- signif(qgamma(pnorm(w3), shape = SH_sim, scale = SC.sim), digits = 4)
            }
            if(verbose && d %% 2 == 0) cat(paste0("\r Realization ", i, ": ", d, "/", ncol(daily_covariates), ". Retries: ", temps_retries, '       '))

            if(d %% 30 == 0) invisible(gc())
        }

        simulated_climate[1, , , 'prcp'][simulated_climate[1, , , 'prcp'] < model$control$prcp_occurrence_threshold] <- 0

        simulated_climate
    }


    attr(gen_climate, 'simulation_coefficients') <- list(
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
