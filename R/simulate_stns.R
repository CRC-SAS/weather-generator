
#' @title Transform simulation result to tibble
#' @description This function transforms the simulation result, a named array, to a tibble.
#' @export
sim_result_to_tibble <- function(sim_result, nsim_to_extract = 1) {

    if (any(class(sim_result) != c("array", "glmwgen.climate")))
        return (NULL)

    dates    <- colnames(sim_result)
    stations <- rownames(sim_result[nsim_to_extract, 1, ,])

    prcp <- tibble::tibble(date = dates) %>%
        dplyr::bind_cols(tibble::as_tibble(sim_result[nsim_to_extract, , , "prcp"])) %>%
        tidyr::gather(key="station", value="prcp", -date)
    tn   <- tibble::tibble(date = dates) %>%
        dplyr::bind_cols(tibble::as_tibble(sim_result[nsim_to_extract, , , "tn"])) %>%
        tidyr::gather(key="station", value="tn", -date)
    tx   <- tibble::tibble(date = dates) %>%
        dplyr::bind_cols(tibble::as_tibble(sim_result[nsim_to_extract, , , "tx"])) %>%
        tidyr::gather(key="station", value="tx", -date)

    result   <- tidyr::crossing(dates, stations) %>%
        dplyr::rename(date = dates, station = stations) %>%
        dplyr::inner_join(prcp, by = c("date","station")) %>%
        dplyr::inner_join(tn, by = c("date","station")) %>%
        dplyr::inner_join(tx, by = c("date","station"))

    return (result)
}


generate_simulation_matrix <- function(daily_covariates, actual_date_index, model_lags,
                                       sim_start_date, start_climatology, simulated_climate) {

    # Las daily_covariates son iguales para todas las estaciones (para una misma fecha), por lo tanto,
    # es necesario repetir la fila para cada una de las estaciones simuladas
    simulation_matrix <- tibble::as_tibble(t(daily_covariates[, actual_date_index])) %>%
        tidyr::crossing(tibble::tibble(station = unique(dplyr::pull(start_climatology, station))))

    # Una vez metidas todas las daily_covariates en simulation_matrix, se agregan datos inciales y previos
    # considerando los lags utilizados en el ajuste, hasta el max lag indicado en model_control
    max_lag <- max(model_lags$tx_lags_to_use, model_lags$tn_lags_to_use, model_lags$prcp_lags_to_use)

    for (i in 1:max_lag) {

        # Al manejar los lags hay dos situaciones; 1- la inicial, cuando los datos previos se toman de
        # start_climatology (dato generado en el ajuste); 2- cuando ya es posible tomar, como previos,
        # datos simulados previamente, en cuyo caso los datos previos se toman de simulated_climate.
        if (actual_date_index <= max_lag) {
            previous_date        <- as.Date(sim_start_date) - i
            previous_climatology <- start_climatology %>%
                filter(month == lubridate::month(previous_date), day == lubridate::day(previous_date))
        } else {
            previous_index       <- actual_date_index - i
            previous_sim_climate <- simulated_climate[1, previous_index, , ]
            if (!is.matrix(previous_sim_climate))
                previous_sim_climate <- t(previous_sim_climate)
            previous_climatology <- tibble::as_tibble(previous_sim_climate)
        }

        # El nombre de las columnas a ser agregadas varía de acuerdo a la cantida de lags
        if (i == 1) sim_matrix_cols <- list(previous_oc = "prcp_occ_prev",
                                            previous_tn = "tn_prev",
                                            previous_tx = "tx_prev",
                                            prcp = "prcp", prcp_occ = "prcp_occ")
        else sim_matrix_cols <- list(previous_oc = paste0("prcp_occ_prev.minus.",i),
                                     previous_tn = paste0("tn_prev.minus.",i),
                                     previous_tx = paste0("tx_prev.minus.",i))

        # Finalmente, aquí se carga simulation_matrix con los datos apropiados
        if (i <= model_lags$prcp_lags_to_use) {
            simulation_matrix    <- simulation_matrix %>%
                dplyr::mutate(!!sim_matrix_cols$previous_oc := as.integer(dplyr::pull(previous_climatology, 'prcp') > 0))
        }
        if (i <= model_lags$tn_lags_to_use) {
            simulation_matrix    <- simulation_matrix %>%
                dplyr::mutate(!!sim_matrix_cols$previous_tn := dplyr::pull(previous_climatology, 'tn'))
        }
        if (i <= model_lags$tx_lags_to_use) {
            simulation_matrix    <- simulation_matrix %>%
                dplyr::mutate(!!sim_matrix_cols$previous_tx := dplyr::pull(previous_climatology, 'tx'))
        }
        if (i == 1) {
            # Estos son datos del día, pero simulados, por lo tanto, en este punto de la ejecución aún
            # no están disponibles. Van a ser calculados más adelante en la ejecución de la simulación!
            simulation_matrix    <- simulation_matrix %>%
                dplyr::mutate(!!sim_matrix_cols$prcp := NA, !!sim_matrix_cols$prcp_occ := NA)
        }

    }

    return(simulation_matrix)

}


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

    if(class(object) != 'glmwgen') stop(paste('Received a model of class', class(object), 'and a model of class "glmwgen" was expected.'))

    simulation_locations <- model$stations

    if(!is.null(seed)) set.seed(seed)

    if(!('proj4string' %in% names(attributes(simulation_locations)))) {
        warning('simulation_locations is not a spacial points object, attempting to convert it with stations projection string.')
        simulation_locations <- SpatialPoints(simulation_locations, proj4string=model$stations_proj4string)
    }

    if(end_date <= start_date) stop('End date should be greater than start date')

    if(nsim < 1) stop('Number of simulations should be greater than one')

    if(identical(control$random_fields_method, glmwgen:::random_field_noise)) {
        warning("There isn't possible to use random_field_noise function to calculate noises when simulated points are the coordinates of stations.\n")
    }

    if(!identical(control$random_fields_method, glmwgen:::rnorm_noise) && (nrow(simulation_locations) == 1)) {
        warning("Only the function rnorm_noise can be used to calculate noises for simulate only one station!\n")
        control$random_fields_method <- glmwgen:::rnorm_noise
    }

    if(!identical(control$random_fields_method, glmwgen:::mvrnorm_noise) && (nrow(simulation_locations) > 1)) {
        warning("Only the function mvnorm_noise can be used to calculate noises for simulate multiple stations!\n")
        control$random_fields_method <- glmwgen:::mvrnorm_noise
    }

    #########################################

    simulation_dates <- data.frame(date = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "days"))

    if (is.null(control$Rt)) {
        # TODO: shouldn't this be = 1 once we start simulating?
        Rt <- seq(from = -1, to = 1, length.out = nrow(simulation_dates))
        # Rt <- 0
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
               six_months_fraction = year_fraction / 2, ct.six = cos(six_months_fraction), st.six = sin(six_months_fraction),
               three_months_fraction = year_fraction / 4, ct.three = cos(three_months_fraction), st.three = sin(three_months_fraction),
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

    daily_covariates <- as.matrix(t(simulation_dates[, !colnames(simulation_dates) %in% c('date', 'year', 'year_fraction', 'six_months_fraction',
                                                                                          'three_months_fraction', 'month', 'season')]))

    daily_covariates <- rbind(1, daily_covariates)
    rownames(daily_covariates)[1] <- "(Intercept)"


    simulation_dates <- simulation_dates %>% dplyr::select(date, year, month)

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
    rownames(coefocc_sim) <- rownames(model$coefficients$coefocc)
    colnames(coefmin_sim) <- colnames(model$coefficients$coefmin)
    rownames(coefmin_sim) <- rownames(model$coefficients$coefmin)
    colnames(coefmax_sim) <- colnames(model$coefficients$coefmax)
    rownames(coefmax_sim) <- rownames(model$coefficients$coefmax)

    SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(sp::coordinates(model$stations))))
    colnames(SC) <- model$stations$id

    for (station in as.character(model$stations$id)) {
        gamma_coef <- model$gamma[[station]]$coef
        SC[, station] <- exp(apply(daily_covariates[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T)) / SH[station]
    }

    #########################################

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

    realizations_seeds <- ceiling(runif(min = 1, max = 10000000, n = nsim))

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

        model_lags <- list(prcp_lags_to_use = model$control$prcp_lags_to_use,
                           tn_lags_to_use = model$control$tn_lags_to_use,
                           tx_lags_to_use = model$control$tx_lags_to_use)

        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {

            simulation_matrix <- glmwgen:::generate_simulation_matrix(daily_covariates, d, model_lags, start_date,
                                                                      model$start_climatology, simulated_climate)

            month_number <- simulation_dates$month[d]

            # Se simula la ocurrencia de precipitación (prcp_occ)
            prcp_occ_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_occ = rowSums(simulation_matrix[, colnames(coefocc_sim)] * coefocc_sim),
                occ_noise = noise_generator$generate_noise(month_number, 'prcp'),
                prcp_occ = as.integer((mu_occ + occ_noise) > 0)
            )

            # Se guarda prcp_occ en simulation_matrix (porque se usa para simular temperatura)
            simulation_matrix$prcp_occ <- dplyr::pull(prcp_occ_sim, prcp_occ)

            # Se calcula una cantidad de mm para las estaciones con lluvia (prcp_amt)
            prcp_amt_sim <- prcp_occ_sim %>%
                dplyr::select(station, prcp_occ) %>%
                dplyr::mutate(
                    amt_noise = noise_generator$generate_noise(month_number, 'prcp'),
                    amt_value = signif(qgamma(pnorm(amt_noise), shape = SH, scale = SC[d, ]), digits = 4)
                ) %>%
                dplyr::mutate(is_valid = prcp_occ == 1 && amt_value >= model$control$prcp_occurrence_threshold) %>%
                dplyr::mutate(prcp_amt = ifelse(is_valid, amt_value, 0))

            # Se guarda prcp (cantidad de mm) en simulation_matrix (porque se usa para simular temperatura)
            simulation_matrix$prcp <- dplyr::pull(prcp_amt_sim, prcp_amt)

            # Se simula la temperatura mínima (tn_sim)
            tn_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_tn = rowSums(simulation_matrix[, colnames(coefmin_sim)] * coefmin_sim),
                tn_noise = noise_generator$generate_noise(month_number, 'tn'),
                tn = signif(mu_tn + tn_noise, digits = 4)
            )

            # Se simula la temperatura máxima (tx_sim)
            tx_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_tx = rowSums(simulation_matrix[, colnames(coefmax_sim)] * coefmax_sim),
                tx_noise = noise_generator$generate_noise(month_number, 'tx'),
                tx = signif(mu_tx + tx_noise, digits = 4)
            )

            # Se verifica que tx y tn sean válidos (sino son válidos se los recalcula)
            daily_retries <- 0
            while (min(tx_sim$tx - tn_sim$tn) < control$minimum_temperatures_difference_threshold && daily_retries < 100) {
                temps_retries <- temps_retries + 1
                daily_retries <- daily_retries + 1
                tn_sim <- tn_sim %>% dplyr::mutate(
                    tn_noise = noise_generator$generate_noise(month_number, 'tn'),
                    tn = signif(mu_tn + tn_noise, digits = 4))
                tx_sim <- tx_sim %>% dplyr::mutate(
                    tx_noise = noise_generator$generate_noise(month_number, 'tx'),
                    tx = signif(mu_tx + tx_noise, digits = 4))
            }
            if(daily_retries >= 100) {
                if(verbose) cat('Failed to simulate random noise that doesn\'t violate the constraint of max. temp. > min. temp.')
                return(NULL)
            }

            # Se guardan prcp, tn y tx en el array de salida (el resultado de la simulación)
            simulated_climate[1, d, , 'prcp'] <- dplyr::pull(prcp_amt_sim, prcp_amt)
            simulated_climate[1, d, , 'tn']   <- dplyr::pull(tn_sim, tn)
            simulated_climate[1, d, , 'tx']   <- dplyr::pull(tx_sim, tx)

            # Se reporta el estado de la simulación
            if(verbose && d %% 2 == 0) cat(paste0("\r Realization ", i, ": ", d, "/", ncol(daily_covariates), ". Retries: ", temps_retries, '       '))
            if(d %% 30 == 0) invisible(gc())
        }

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
