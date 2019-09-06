

generate_stns_simulation_matrix <- function(daily_covariates, actual_date, actual_date_index,
                                            sim_start_date, start_climatology, simulated_climate) {

    simulation_matrix <- daily_covariates %>% dplyr::filter(date == actual_date)

    # Una vez metidas todas las daily_covariates en simulation_matrix, se agregan datos inciales y previos.
    # Al hacerlo, se presentan dos situaciones posibles;
    if (actual_date_index == 1) {
        # 1- la inicial, cuando los datos previos se toman de start_climatology (datos generados en el ajuste);
        previous_date        <- as.Date(sim_start_date) - 1
        previous_climatology <- start_climatology %>%
            filter(month == lubridate::month(previous_date), day == lubridate::day(previous_date))
    } else {
        # 2- cuando ya es posible tomar, como previos, datos simulados en pasos anteriores,
        # en cuyo caso los datos previos se toman de simulated_climate.
        previous_index       <- actual_date_index - 1
        previous_sim_climate <- simulated_climate[1, previous_index, , ]
        if (!is.matrix(previous_sim_climate))
            previous_sim_climate <- t(previous_sim_climate)
        previous_climatology <- tibble::as_tibble(previous_sim_climate)
    }

    # Se definen las nuevas columnas a ser agregadas a simulation_matrix
    sim_matrix_cols <- list(previous_oc = "prcp_occ_prev",
                            previous_tn = "tn_prev", tn = "tn",
                            previous_tx = "tx_prev", tx = "tx",
                            prcp = "prcp", prcp_occ = "prcp_occ")

    # Estos son datos de días previos, de start_climatology (datos generados en el ajuste) o de
    # simulated_climate (datos generados en pasos previos), según corresponda.
    simulation_matrix    <- simulation_matrix %>%
        dplyr::mutate(!!sim_matrix_cols$previous_oc := as.integer(dplyr::pull(previous_climatology, 'prcp') > 0))
    simulation_matrix    <- simulation_matrix %>%
        dplyr::mutate(!!sim_matrix_cols$previous_tn := dplyr::pull(previous_climatology, 'tn'))
    simulation_matrix    <- simulation_matrix %>%
        dplyr::mutate(!!sim_matrix_cols$previous_tx := dplyr::pull(previous_climatology, 'tx'))

    # Estos son datos del día, pero simulados, por lo tanto, en este punto de la ejecución aún
    # no están disponibles. Van a ser calculados más adelante en la ejecución de la simulación!
    simulation_matrix    <- simulation_matrix %>%
        dplyr::mutate(!!sim_matrix_cols$tn := NA, !!sim_matrix_cols$tx := NA) %>%
        dplyr::mutate(!!sim_matrix_cols$prcp := NA, !!sim_matrix_cols$prcp_occ := NA)

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


    #########################################


    unique_stations <- unique(model$stations$id) # en lugar de matching_stations

    # Save original coordinates.
    simulation_coordinates <- sp::coordinates(simulation_locations)
    rownames(simulation_coordinates) <- model$stations$id


    #########################################


    simulation_dates <- data.frame(date = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "days")) %>%
        dplyr::mutate(year   = lubridate::year(date),
                      month  = lubridate::month(date),
                      day    = lubridate::day(date),
                      doy    = lubridate::yday(date),
                      season = lubridate::quarter(date, fiscal_start = 12))


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
                      Rt = sort(rep(Rt, times = length(simulation_locations))),
                      row_num = row_number()) %>%
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
        dplyr::rename(sum_prcp = st_covariates,
                      mean_tx  = smx_covariates,
                      mean_tn  = smn_covariates) %>%
        dplyr::left_join(do.call("rbind", model$residuals), by = c("station", "date")) %>%
        dplyr::select(station, date, dplyr::everything(), -month, -day, -dplyr::ends_with("fraction"))


    estadisticos_umbrales <- tibble::tibble(station = integer(), date = lubridate::ymd(), prcp_occ = integer(),
                                            max.range = integer(), min.range = integer())
    for ( station_id in names(model$estadisticos.umbrales) ) {
        for ( mes in unique(simulation_dates$month) ) {
            for ( prcp in c("con.prcp", "sin.prcp") ) {
                estadisticos_umbrales <- estadisticos_umbrales %>%
                    dplyr::bind_rows(tibble::tibble(
                        station   = as.integer(station_id),
                        date      = simulation_dates %>% dplyr::filter(month == mes) %>% dplyr::pull(date),
                        prcp_occ  = ifelse(prcp == "con.prcp", as.integer(T), as.integer(F)),
                        max.range = model$estadisticos.umbrales[[station_id]][[mes]][[prcp]][["max.range"]],
                        min.range = model$estadisticos.umbrales[[station_id]][[mes]][[prcp]][["min.range"]]))
            }
        }
    }


    #########################################


    # Register a sequential backend if the user didn't register a parallel
    # in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        foreach::registerDoSEQ()
    }

    combine_function <- function(...) {
        invisible(gc())
        abind::abind(..., along = 1)
    }

    global_retries <- 0
    realizations_seeds <- ceiling(runif(min = 1, max = 10000000, n = nsim))

    # gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %dopar% {
    gen_climate <- foreach(i = 1:nsim, .combine = combine_function, .multicombine = control$multicombine, .packages = c('sp')) %dopar% {

        simulated_climate <- array(data = 0.0, dim = c(1, nrow(simulation_dates), nrow(simulation_coordinates), 3))
        dimnames(simulated_climate)[2] <- list('dates' = format(simulation_dates$date, '%Y-%m-%d'))
        dimnames(simulated_climate)[3] <- list('coordinates' = rownames(simulation_coordinates))
        dimnames(simulated_climate)[4] <- list('variables' = c('tx', 'tn', 'prcp'))

        ruidos_aleatorios <- tibble::tibble(station = integer(), date = lubridate::ymd(), prcp_occ = integer(),
                                            tx.noise = double(), tn.noise = double())
        for ( station_id in names(model$estadisticos.residuos) ) {
            for ( mes in unique(simulation_dates$month) ) {
                for ( prcp in c("con.prcp", "sin.prcp") ) {
                    set.seed(realizations_seeds[i]) # para cuando necesitamos repetir resultados
                    ruidos.tx.tn <- MASS::mvrnorm(n = nrow(simulation_dates %>% dplyr::filter(month == mes)),
                                                  mu = model$estadisticos.residuos[[station_id]][[mes]][[prcp]][["media"]],
                                                  Sigma = model$estadisticos.residuos[[station_id]][[mes]][[prcp]][["matriz.covarianza"]],
                                                  empirical = TRUE) %>% magrittr::set_colnames(c("tx.noise", "tn.noise"))
                    ruidos_aleatorios <- ruidos_aleatorios %>%
                        dplyr::bind_rows(tibble::tibble(
                            station   = as.integer(station_id),
                            date      = simulation_dates %>% dplyr::filter(month == mes) %>% dplyr::pull(date),
                            prcp_occ  = ifelse(prcp == "con.prcp", as.integer(T), as.integer(F)),
                            tx.noise  = ruidos.tx.tn[,"tx.noise"],
                            tn.noise  = ruidos.tx.tn[,"tn.noise"]))
                }
            }
        }

        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {
            actual_date <- simulation_dates$date[d]

            simulation_matrix <- glmwgen:::generate_stns_simulation_matrix(daily_covariates, actual_date, d, start_date,
                                                                           model$start_climatology, simulated_climate)

            month_number <- simulation_dates$month[d]

            # Se simula la ocurrencia de precipitación (prcp_occ)
            prcp_occ_sim <- foreach::foreach(station_id = unique_stations, .multicombine = T, .combine = dplyr::bind_rows) %dopar% {
                set.seed(realizations_seeds[i]) # para cuando necesitamos repetir resultados
                #
                prcp_occ_fit <- model$gam_fits[[as.character(station_id)]]$prcp_occ_fit
                stn_sim_mtrx <- simulation_matrix %>% dplyr::filter(station == as.integer(station_id))
                #
                result_predict <- mgcv::predict.gam(object = prcp_occ_fit, newdata = stn_sim_mtrx)
                result_rnorm   <- stats::rnorm(n = 1, mean = 0, sd = 1)
                # occurrence is mgcv::predict.gam(prcp_occ_fit) + rnorm(mean=0, sd=1)
                tibble::tibble(station = station_id, date = actual_date,
                               prcp_occ_simulated = as.integer(result_predict + result_rnorm > 0))
            }

            # Se guarda prcp_occ en simulation_matrix (porque se usa para simular temperaturas)
            simulation_matrix <- simulation_matrix %>%
                dplyr::left_join(prcp_occ_sim, by = c("station", "date"), suffix = c("", "_simulated")) %>%
                dplyr::mutate(prcp_occ = prcp_occ_simulated) %>%
                dplyr::select(-prcp_occ_simulated)


            # Se calcula una cantidad de mm para las estaciones con lluvia (prcp_amt)
            prcp_amt_sim <- foreach::foreach(station_id = unique_stations, .multicombine = T, .combine = dplyr::bind_rows) %dopar% {
                set.seed(realizations_seeds[i]) # para cuando necesitamos repetir resultados
                #
                if (prcp_occ_sim %>% dplyr::filter(station == as.integer(station_id)) %>% dplyr::pull(prcp_occ_simulated) %>% magrittr::not())
                    return (tibble::tibble(station = station_id, date = actual_date, prcp_amt_simulated = 0))
                #
                prcp_amt_fit <- model$gam_fits[[as.character(station_id)]]$prcp_amt_fit[[month_number]]
                stn_sim_mtrx <- simulation_matrix %>% dplyr::filter(station == as.integer(station_id))
                #
                result_predict <- mgcv::predict.gam(object = prcp_amt_fit, newdata = stn_sim_mtrx)
                result_pnorm   <- stats::pnorm(stats::rnorm(n = 1, mean = 0, sd = 1))
                # Estimación del parametro de forma (alpha_amt) y el de escala (beta_amt)
                alpha_amt <- MASS::gamma.shape(prcp_amt_fit)$alpha
                beta_amt  <- exp(result_predict) / alpha_amt
                #
                tibble::tibble(station = station_id, date = actual_date,
                               prcp_amt_simulated = stats::qgamma(result_pnorm, shape = alpha_amt, scale = beta_amt))
            }

            # Se guarda prcp (cantidad de mm) en simulation_matrix (porque se usa para simular temperaturas)
            simulation_matrix <- simulation_matrix %>%
                dplyr::left_join(prcp_amt_sim, by = c("station", "date"), suffix = c("", "_simulated")) %>%
                dplyr::mutate(prcp_amt = prcp_amt_simulated) %>%
                dplyr::select(-prcp_amt_simulated)


            # Se simula la temperatura mínima (tn_sim)
            tn_sim <- foreach::foreach(station_id = unique_stations, .multicombine = T, .combine = dplyr::bind_rows) %dopar% {
                set.seed(realizations_seeds[i]) # para cuando necesitamos repetir resultados
                #
                tn_sim_fit         <- model$gam_fits[[as.character(station_id)]]$tn_fit
                stn_sim_mtrx       <- simulation_matrix %>% dplyr::filter(station == as.integer(station_id))
                prcp_occ_simulated <- prcp_occ_sim %>% dplyr::filter(station == as.integer(station_id)) %>% dplyr::pull(prcp_occ_simulated)
                #
                result_predict  <- mgcv::predict.gam(object = tn_sim_fit, newdata = stn_sim_mtrx)
                ruido_aleatorio <- ruidos_aleatorios %>%
                    dplyr::filter(station == as.integer(station_id), date == actual_date, prcp_occ == prcp_occ_simulated) %>%
                    dplyr::pull(tn.noise)
                #
                tibble::tibble(station = station_id, date = actual_date, tn_simulated = result_predict + ruido_aleatorio)
            }

            # Se guarda tn en simulation_matrix (porque se usa para simular tx)
            simulation_matrix <- simulation_matrix %>%
                dplyr::left_join(tn_sim, by = c("station", "date"), suffix = c("", "_simulated")) %>%
                dplyr::mutate(tn_sim = tn_simulated) %>%
                dplyr::select(-tn_simulated)


            # Se simula la temperatura máxima (tx_sim)
            tx_sim <- foreach::foreach(station_id = unique_stations, .multicombine = T, .combine = dplyr::bind_rows) %dopar% {
                set.seed(realizations_seeds[i]) # para cuando necesitamos repetir resultados
                #
                tx_sim_fit         <- model$gam_fits[[as.character(station_id)]]$tx_fit
                stn_sim_mtrx       <- simulation_matrix %>% dplyr::filter(station == as.integer(station_id))
                prcp_occ_simulated <- prcp_occ_sim %>% dplyr::filter(station == as.integer(station_id)) %>% dplyr::pull(prcp_occ_simulated)
                #
                result_predict  <- mgcv::predict.gam(object = tx_sim_fit, newdata = stn_sim_mtrx)
                ruido_aleatorio <- ruidos_aleatorios %>%
                    dplyr::filter(station == as.integer(station_id), date == actual_date, prcp_occ == prcp_occ_simulated) %>%
                    dplyr::pull(tx.noise)
                #
                tibble::tibble(station = station_id, date = actual_date, tx_simulated = result_predict + ruido_aleatorio)
            }

            # Se guarda tx en simulation_matrix (solo para mejorar los debugs, no se usa en ningún cálculo posterior)
            simulation_matrix <- simulation_matrix %>%
                dplyr::left_join(tx_sim, by = c("station", "date"), suffix = c("", "_simulated")) %>%
                dplyr::mutate(tx_sim = tx_simulated) %>%
                dplyr::select(-tx_simulated)


            # Se establecen los parametros de control de las temperaturas generadas
            t_ctrl <- tidyr::crossing(station = unique_stations, date = actual_date) %>%
                dplyr::inner_join(prcp_occ_sim, by = c("station", "date")) %>%
                dplyr::inner_join(tx_sim,       by = c("station", "date")) %>%
                dplyr::inner_join(tn_sim,       by = c("station", "date")) %>%
                dplyr::left_join(estadisticos_umbrales, by = c("station", "date", "prcp_occ_simulated" = "prcp_occ")) %>%
                dplyr::mutate(te = tx_simulated - tn_simulated) %>%
                dplyr::select(station, date, tx = tx_simulated, tn = tn_simulated, te_min = min.range, te, te_max = max.range)


            # Se verifica que tx y tn sean válidos (sino son válidos se los recalcula)
            daily_retries <- 0
            while ( daily_retries < 100 && (any(t_ctrl$tx < t_ctrl$tn) || any(t_ctrl$te > t_ctrl$te_max) || any(t_ctrl$te < t_ctrl$te_min)) ) {
                temps_retries  <- temps_retries  + 1
                daily_retries  <- daily_retries  + 1
                global_retries <- global_retries + 1

                stns_a_recalc <- t_ctrl %>% dplyr::filter(tx < tn | te > te_max | te < te_min) %>% dplyr::pull(station)

                if(verbose) cat("\n"); cat("OLD temp_control"); print(t_ctrl);

                for (station_id in stns_a_recalc) {
                    tx_sim_fit <- model$gam_fits[[as.character(station_id)]]$tx_fit
                    tn_sim_fit <- model$gam_fits[[as.character(station_id)]]$tn_fit

                    stn_sim_mtrx <- simulation_matrix %>% dplyr::filter(station == as.integer(station_id))

                    # Generar nuevas temperaturas
                    set.seed(realizations_seeds[i] + daily_retries) # para cuando necesitamos repetir resultados
                    new_tx <- mgcv::predict.gam(tx_sim_fit, newdata = stn_sim_mtrx) + rnorm(n = 1, mean = 0, sd = stn_sim_mtrx$sd.tx_residuals)
                    set.seed(realizations_seeds[i] + daily_retries) # para cuando necesitamos repetir resultados
                    new_tn <- mgcv::predict.gam(tn_sim_fit, newdata = stn_sim_mtrx) + rnorm(n = 1, mean = 0, sd = stn_sim_mtrx$sd.tn_residuals)

                    # Actualizar temperaturas simuladas
                    tx_sim <- tx_sim %>% dplyr::mutate(tx_simulated = ifelse(station == as.integer(station_id), new_tx, tx_simulated))
                    tn_sim <- tn_sim %>% dplyr::mutate(tn_simulated = ifelse(station == as.integer(station_id), new_tn, tn_simulated))

                    # Actualizar los parametros de control de las temperaturas generadas
                    t_ctrl <- t_ctrl %>% dplyr::mutate(tx = ifelse(station == as.integer(station_id), new_tx, tx),
                                                       tn = ifelse(station == as.integer(station_id), new_tn, tn),
                                                       te = ifelse(station == as.integer(station_id), tx - tn, te))
                }
                if(verbose) cat("NEW temp_control"); print(t_ctrl);
            }
            if(daily_retries >= 100) {
                if(verbose) cat('Failed to simulate random noise that doesn\'t violate the constraint of max. temp. > min. temp.')
                return(NULL)
            }

            # Se guardan prcp, tn y tx en el array de salida (el resultado de la simulación)
            simulated_climate[1, d, , 'prcp'] <- dplyr::pull(prcp_amt_sim, prcp_amt_simulated)
            simulated_climate[1, d, , 'tn']   <- dplyr::pull(tn_sim,       tn_simulated)
            simulated_climate[1, d, , 'tx']   <- dplyr::pull(tx_sim,       tx_simulated)

            # Se reporta el estado de la simulación
            if(verbose && d %% 2 == 0) cat(paste0("\r Realization ", i, ": ", d, "/", nrow(simulation_dates), ". Retries: ", temps_retries, '       '))
            if(d %% 30 == 0) invisible(gc())
        }

        simulated_climate
    }

    if(verbose) cat("\n")  # Para que los warnings y errores se impriman debajo de Realization ..., y no al lado!


    attr(gen_climate, 'seasonal_covariates') <- list(
        'tx' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, smx_covariates),  # smx_covariates,
        'tn' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, smn_covariates),  # smn_covariates,
        'prcp' = seasonal_covariates %>% dplyr::filter(year %in% years_in_sim_dates) %>%
            dplyr::select(station, year, season, st_covariates)    # st_covariates
    )

    attr(gen_climate, 'tx_tn_retries') <- global_retries
    attr(gen_climate, 'realizations_seeds') <- realizations_seeds
    attr(gen_climate, 'simulation_coordinates') <- simulation_coordinates

    attr(gen_climate, 'model') <- model

    class(gen_climate) <- c(class(gen_climate), 'glmwgen.climate')

    gen_climate
}
