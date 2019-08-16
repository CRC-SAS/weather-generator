

generate_stns_simulation_matrix <- function(daily_covariates, actual_date, actual_date_index, model_lags,
                                            sim_start_date, start_climatology, simulated_climate) {

    simulation_matrix <- daily_covariates %>% dplyr::filter(date == actual_date)

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
                                            previous_tn = "tn_prev", tn = "tn",
                                            previous_tx = "tx_prev", tx = "tx",
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
                dplyr::mutate(!!sim_matrix_cols$tn := NA, !!sim_matrix_cols$tx := NA) %>%
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
                             control = glmwgen:::glmwgen_simulation_control(), verbose = T) {
    model <- object

    if(class(object) != 'glmwgen') stop(paste('Received a model of class', class(object), 'and a model of class "glmwgen" was expected.'))

    simulation_locations <- model$stations

    if(!is.null(seed)) set.seed(seed)

    if(!('proj4string' %in% names(attributes(simulation_locations)))) {
        warning('simulation_locations is not a spacial points object, attempting to convert it with stations projection string.')
        simulation_locations <- SpatialPoints(simulation_locations, proj4string=model$stations_proj4string)
    }

    if(end_date <= start_date) stop('End date should be greater than start date')

    years_in_sim_dates <- base::seq.int(lubridate::year(start_date), lubridate::year(end_date))
    years_in_senal_cov <- dplyr::distinct(model$seasonal, year) %>% dplyr::pull()
    if (!all(is.element(years_in_sim_dates, years_in_senal_cov)))
        stop("Simulation years aren't in model$seasonal!")

    if(nsim < 1) stop('Number of simulations should be greater than one')

    if(identical(control$random_fields_method, glmwgen:::random_field_noise)) {
        warning("There isn't possible to use random_field_noise function to calculate noises when simulated points are the coordinates of stations.\n")
    }

    if(!identical(control$random_fields_method, glmwgen:::rnorm_noise) && (nrow(simulation_locations) == 1)) {
        warning("Only the function rnorm_noise can be used to calculate noises for simulate only one station! Switching to it!\n")
        control$random_fields_method <- glmwgen:::rnorm_noise
    }

    if(!identical(control$random_fields_method, glmwgen:::mvrnorm_noise) && (nrow(simulation_locations) > 1)) {
        warning("Only the function mvnorm_noise can be used to calculate noises for simulate multiple stations! Switching to it!\n")
        control$random_fields_method <- glmwgen:::mvrnorm_noise
    }

    #########################################

    simulation_dates <- data.frame(date = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "days")) %>%
        dplyr::mutate(year = lubridate::year(date), month = lubridate::month(date), day = lubridate::day(date))

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

    seasonal_covariates <- model$seasonal  %>%
        dplyr::rename(st_covariates  = dplyr::ends_with("prcp"),
                      smx_covariates = dplyr::ends_with("tx"),
                      smn_covariates = dplyr::ends_with("tn"))

    daily_covariates    <- simulation_dates %>%
        tidyr::crossing(station = simulation_locations$id) %>%
        dplyr::mutate("(Intercept)" = 1,
                      year = lubridate::year(date),
                      year_fraction = 2 * pi * lubridate::yday(date) / ifelse(lubridate::leap_year(date), 366, 365),
                      ct = cos(year_fraction),
                      st = sin(year_fraction),
                      six_months_fraction = year_fraction / 2,
                      ct.six = cos(six_months_fraction),
                      st.six = sin(six_months_fraction),
                      Rt = sort(rep(Rt, times = length(simulation_locations))),
                      season = ceiling(lubridate::month(date)/3)) %>%
        dplyr::left_join(seasonal_covariates, by = c("station", "year", "season")) %>%
        dplyr::mutate(ST1  = dplyr::if_else(season == 1, st_covariates,  0),
                      ST2  = dplyr::if_else(season == 2, st_covariates,  0),
                      ST3  = dplyr::if_else(season == 3, st_covariates,  0),
                      ST4  = dplyr::if_else(season == 4, st_covariates,  0),
                      SMX1 = dplyr::if_else(season == 1, smx_covariates, 0),
                      SMX2 = dplyr::if_else(season == 2, smx_covariates, 0),
                      SMX3 = dplyr::if_else(season == 3, smx_covariates, 0),
                      SMX4 = dplyr::if_else(season == 4, smx_covariates, 0),
                      SMN1 = dplyr::if_else(season == 1, smn_covariates, 0),
                      SMN2 = dplyr::if_else(season == 2, smn_covariates, 0),
                      SMN3 = dplyr::if_else(season == 3, smn_covariates, 0),
                      SMN4 = dplyr::if_else(season == 4, smn_covariates, 0)) %>%
        dplyr::select(station, date, dplyr::everything(), -month, -day, -dplyr::ends_with("fraction"))


    #########################################

    stations <- seq_len(nrow(simulation_locations))  # en lugar de matching_stations

    # Save original coordinates.
    simulation_coordinates <- sp::coordinates(simulation_locations)
    rownames(simulation_coordinates) <- model$stations$id


    coefocc_sim <- matrix(model$coefficients$coefocc[stations, ], nrow = length(stations), byrow = F)
    coefmin_sim <- matrix(model$coefficients$coefmin[stations, ], nrow = length(stations), byrow = F)
    coefmax_sim <- matrix(model$coefficients$coefmax[stations, ], nrow = length(stations), byrow = F)

    colnames(coefocc_sim) <- colnames(model$coefficients$coefocc)
    rownames(coefocc_sim) <- rownames(model$coefficients$coefocc)
    colnames(coefmin_sim) <- colnames(model$coefficients$coefmin)
    rownames(coefmin_sim) <- rownames(model$coefficients$coefmin)
    colnames(coefmax_sim) <- colnames(model$coefficients$coefmax)
    rownames(coefmax_sim) <- rownames(model$coefficients$coefmax)


    # Gamma shape
    SH <- c()
    for (station in as.character(model$stations$id)) SH[station] <- model$gamma[[station]]$alpha
    SH_sim <- SH[stations]


    SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(sp::coordinates(model$stations))))
    colnames(SC) <- model$stations$id

    for (station_id in as.character(model$stations$id)) {
        gamma_coef <- model$gamma[[station_id]]$coef
        col_names  <- as.character(sort(simulation_dates$date))
        covariates <- daily_covariates %>%
            dplyr::filter(station == station_id) %>% dplyr::arrange(date) %>%
            dplyr::select(-station, -date) %>% base::t() %>% base::as.matrix() %>%
            magrittr::set_colnames(glue::glue("{station_id}.{col_names}"))
        SC[, station_id] <- exp(apply(covariates[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T)) / SH[station_id]
    }

    SC_sim <- SC

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
            actual_date <- simulation_dates$date[d]

            simulation_matrix <- glmwgen:::generate_stns_simulation_matrix(daily_covariates, actual_date, d, model_lags, start_date,
                                                                           model$start_climatology, simulated_climate)

            month_number <- simulation_dates$month[d]

            # Se simula la ocurrencia de precipitación (prcp_occ)
            prcp_occ_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_occ = rowSums(simulation_matrix[, colnames(coefocc_sim)] * coefocc_sim),
                occ_noise = noise_generator$generate_noise(month_number, 'prcp'),
                prcp_occ = as.integer((mu_occ + occ_noise) > 0)
            )

            # Se guarda prcp_occ en simulation_matrix (porque se usa para simular temperaturas)
            simulation_matrix$prcp_occ <- dplyr::pull(prcp_occ_sim, prcp_occ)

            # Se calcula una cantidad de mm para las estaciones con lluvia (prcp_amt)
            prcp_amt_sim <- prcp_occ_sim %>%
                dplyr::select(station, prcp_occ) %>%
                dplyr::mutate(
                    amt_noise = noise_generator$generate_noise(month_number, 'prcp'),
                    amt_value = signif(qgamma(pnorm(amt_noise), shape = SH_sim, scale = SC_sim[d, ]), digits = 4)
                ) %>%
                dplyr::mutate(is_valid = prcp_occ == 1 && amt_value >= model$control$prcp_occurrence_threshold) %>%
                dplyr::mutate(prcp_amt = ifelse(is_valid, amt_value, 0))

            # Se guarda prcp (cantidad de mm) en simulation_matrix (porque se usa para simular temperaturas)
            simulation_matrix$prcp <- dplyr::pull(prcp_amt_sim, prcp_amt)

            # Se simula la temperatura mínima (tn_sim)
            tn_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_tn = rowSums(simulation_matrix[, colnames(coefmin_sim)] * coefmin_sim),
                tn_noise = noise_generator$generate_noise(month_number, 'tn'),
                tn = signif(mu_tn + tn_noise, digits = 4)
            )

            # Se guarda tn en simulation_matrix (porque se usa para simular tx)
            simulation_matrix$tn <- dplyr::pull(tn_sim, tn)

            # Se simula la temperatura máxima (tx_sim)
            tx_sim <- tibble::tibble(
                station = rownames(simulation_coordinates),
                mu_tx = rowSums(simulation_matrix[, colnames(coefmax_sim)] * coefmax_sim),
                tx_noise = noise_generator$generate_noise(month_number, 'tx'),
                tx = signif(mu_tx + tx_noise, digits = 4)
            )

            # Se guarda tx en simulation_matrix (solo para mejorar los debugs, no se usa en ningún cálculo posterior)
            simulation_matrix$tx <- dplyr::pull(tx_sim, tx)

            # Se establecen los parametros de control de las temperaturas generadas
            t_ctrl <- model$month_params[[month_number]]$temp_ampl %>%
                tibble::rownames_to_column(var = "station") %>%
                dplyr::inner_join(dplyr::select(tx_sim, station, tx), by = "station") %>%
                dplyr::inner_join(dplyr::select(tn_sim, station, tn), by = "station") %>%
                dplyr::mutate(te = tx - tn) %>%
                dplyr::select(station, tx, tn, te_min, te, te_max)


            # Se verifica que tx y tn sean válidos (sino son válidos se los recalcula)
            daily_retries <- 0
            while ( daily_retries < 100 && (any(t_ctrl$tx < t_ctrl$tn) || any(t_ctrl$te > t_ctrl$te_max) || any(t_ctrl$te < t_ctrl$te_min)) ) {
                temps_retries <- temps_retries + 1
                daily_retries <- daily_retries + 1

                stns_a_recalc <- t_ctrl %>% dplyr::filter(tx < tn | te > te_max | te < te_min) %>% dplyr::pull(station)

                new_tn_noises <- noise_generator$generate_noise(month_number, 'tn')
                names(new_tn_noises) <- dplyr::pull(t_ctrl, station)  # a veces faltan los names, ej. cuando se usa rnorm_noise
                new_tx_noises <- noise_generator$generate_noise(month_number, 'tx')
                names(new_tx_noises) <- dplyr::pull(t_ctrl, station)  # a veces faltan los names, ej. cuando se usa rnorm_noise

                tn_sim <- tn_sim %>% dplyr::mutate(
                    tn_noise = dplyr::if_else(station %in% stns_a_recalc, new_tn_noises[station], tn_noise),
                    tn = dplyr::if_else(station %in% stns_a_recalc, signif(mu_tn + tn_noise, digits = 4), tn))
                tx_sim <- tx_sim %>% dplyr::mutate(
                    tx_noise = dplyr::if_else(station %in% stns_a_recalc, new_tx_noises[station], tx_noise),
                    tx = dplyr::if_else(station %in% stns_a_recalc, signif(mu_tx + tx_noise, digits = 4), tx))
                t_ctrl <- t_ctrl %>% dplyr::mutate(tx = tx_sim$tx, tn = tn_sim$tn, te = tx - tn)
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
            if(verbose && d %% 2 == 0) cat(paste0("\r Realization ", i, ": ", d, "/", nrow(simulation_dates), ". Retries: ", temps_retries, '       '))
            if(d %% 30 == 0) invisible(gc())
        }

        simulated_climate
    }

    if(verbose) cat("\n")  # Para que los warnings y errores se impriman debajo de Realization ..., y no al lado!


    attr(gen_climate, 'simulation_coefficients') <- list(
        'tx' = coefmax_sim,
        'tn' = coefmin_sim,
        'occ' = coefocc_sim
    )

    attr(gen_climate, 'seasonal_covariates') <- list(
        'tx' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, smx_covariates),  # smx_covariates,
        'tn' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, smn_covariates),  # smn_covariates,
        'prcp' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, st_covariates)    # st_covariates
    )

    attr(gen_climate, 'realizations_seeds') <- realizations_seeds
    attr(gen_climate, 'simulation_coordinates') <- simulation_coordinates

    attr(gen_climate, 'model') <- model

    class(gen_climate) <- c(class(gen_climate), 'glmwgen.climate')

    gen_climate
}
