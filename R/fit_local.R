
#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
local_fit_control <- function(prcp_occurrence_threshold = 0.1,
                              use_external_seasonal_climate = T,
                              climate_missing_threshold = 0.2,
                              use_covariates = F, avbl_cores = 2,
                              planar_crs_in_metric_coords = 22185) {

    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                use_external_seasonal_climate = use_external_seasonal_climate,
                climate_missing_threshold = climate_missing_threshold,
                use_covariates = use_covariates, avbl_cores = avbl_cores,
                planar_crs_in_metric_coords = planar_crs_in_metric_coords))
}


#' @title Weather model calibration
#' @description Fits a weather model from historical data.
#' @import foreach
#' @import dplyr
#' @export
local_calibrate <- function(climate, stations, seasonal_climate = NULL,
                            control = gamwgen:::local_fit_control(),
                            verbose = F) {

    ## Se crea el objeto a ser retornado al culminar el ajuste!
    model <- list()

    ###############################################################

    if (control$use_external_seasonal_climate & is.null(seasonal_climate))
        stop("If use_external_seasonal_climate is True, seasonal_climate can't be NULL!!")

    if (!control$use_external_seasonal_climate & !is.null(seasonal_climate)) {
        warning('Entered seasonal_climate will be descarted because use_external_seasonal_climate is set as False!')
        seasonal_climate <- NULL
    }

    # Se controlan que los datos recibidos tengan el formato correcto
    gamwgen:::check.fit.input.data(climate, stations, seasonal_climate)

    ###############################################################

    climate <- climate %>%
        dplyr::arrange(station_id, date) %>%
        dplyr::mutate(prcp_occ_threshold = control$prcp_occurrence_threshold) %>%
        dplyr::mutate(prcp_occ = as.integer(prcp >= control$prcp_occurrence_threshold),  # prcp occurrence
                      tipo_dia = factor(ifelse(as.logical(prcp_occ), 'Lluvioso', 'Seco'), levels = c('Lluvioso', 'Seco')),
                      prcp_amt = ifelse(as.logical(prcp_occ), prcp, NA_real_))  # prcp amount/intensity

    invalid_records <- c(which(climate$tmax < climate$tmin), which(climate$prcp < 0))

    if(length(invalid_records) > 0) {
        str_stns <- paste(unique(climate[invalid_records, 'station_id']), collapse = ", " )
        warning(glue::glue('Setting to NA {length(invalid_records)} records with maximum ',
                           'temperature < minimum temperature and/or negative rainfalls. ',
                           'The afected stations are: {str_stns}.'))
        climate[invalid_records, c('tmax', 'tmin', 'prcp')] <- NA
    }

    ###############################################################

    stations <- stations %>%
        sf::st_transform(control$planar_crs_in_metric_coords) %>%
        dplyr::mutate(longitude = sf::st_coordinates(geometry)[,'X'],
                      latitude  = sf::st_coordinates(geometry)[,'Y']) %>%
        dplyr::arrange(station_id)

    # Se obtienen los estaciones en el df con datos climaticos
    unique_stations <- sort(unique(climate$station_id))

    # Se verifique que las estaciones en el df con datos climaticos sean las mismas que en el df de estaciones
    if (!all(unique_stations == stations$station_id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    ###############################################################

    if (control$use_external_seasonal_climate)
        seasonal_climate <- seasonal_climate %>% dplyr::arrange(station_id, year, season)

    t <- proc.time()
    # Si no se recibió seasonal_climate como parámetro, entonces se la calcula internamente!!
    if (is.null(seasonal_climate) | !control$use_external_seasonal_climate)
        seasonal_climate <- gamwgen:::summarise_seasonal_climate(climate, control$climate_missing_threshold)
    tiempo.seasonal_climate <- proc.time() - t

    # Se verifica que hayan covariables suficientes para cubrir todas las fechas en climate,
    # pero el control solo se hace si se van a utilizar covariables en el ajuste!!
    climate_control <- climate %>% dplyr::distinct(station_id, year = lubridate::year(date),
                                                   season = lubridate::quarter(date, fiscal_start = 12))
    seasnl_cli_ctrl <- seasonal_climate %>% dplyr::distinct(station_id, year, season)
    if (!dplyr::all_equal(climate_control, seasnl_cli_ctrl) & control$use_covariates)
        stop("Years in climate and years in seasonal_climate don't match!")

    ###############################################################



    ######################################
    ## AQUI EMPIEZA REALMENTE EL AJUSTE ##
    ######################################


    #################################
    ## Control de tiempo de ejecución
    t.m <- proc.time()


    ###########################################
    ## Se guardan en model los datos de entrada

    model[["control"]] <- control

    model[["stations"]] <- stations

    model[["climate"]] <- climate

    model[["seasonal_data"]] <- seasonal_climate

    model[["crs_used_to_fit"]] <- sf::st_crs(control$planar_crs_in_metric_coords)



    ####################################
    ## Parallelization initialization ##
    ####################################


    ## Variable that indicate if it's necessary to remove
    ## the parallelization configuration
    remove_parallelization_conf <- F

    ## Register a sequential or a parallel backend for the foreach package
    if (control$avbl_cores == 1) {
        foreach::registerDoSEQ()
    } else {
        remove_parallelization_conf <- T
        # Create cluster
        cluster  <- snow::makeCluster(type = "SOCK",
                                      spec = rep('localhost', length.out = control$avbl_cores),
                                      outfile = ifelse(verbose, "", snow::defaultClusterOptions$outfile))
        # Register cluster as backend for the %dopar% function
        doSNOW::registerDoSNOW(cluster)
    }

    ## Register the number of workers to decide
    ## how manage progress bar in child processes
    nworkers <- foreach::getDoParWorkers()



    ###################################
    ## PREPARACIÓN BARRA DE PROGRESO ##
    ###################################


    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = paste0(ifelse(nworkers == 1, " fitting: :what", " fitted stations:"),
                        ifelse(nworkers == 1, " | station: :stn (:c / ", " :c / "),
                        length(unique_stations), ifelse(nworkers == 1, ")", " | :what"),
                        " | progress: :bar :percent (in :elapsed) | eta: :eta"),
        total = (length(unique_stations)*90)+11, clear = FALSE, width= 90, show_after = 0)

    ## For manage the progress bar in parallel executions
    progress_pb <- function(c) {
        pb$tick(90, tokens = list(c = c, what = "fitting" ))
    }



    ############################
    ## PREPARACIÓN DEL AJUSTE ##
    ############################


    ##############
    # Progress Bar
    pb$tick(0, tokens = list(c = 0, what = "init fit", stn = "--"))


    #####################################
    # Control de tiempo de la preparación
    t.p <- proc.time()


    # Group by stations and conf cluster
    models_data <- climate %>%
        dplyr::group_by(station_id)

    # Crear matriz con todas las estaciones
    t <- proc.time()
    models_data <- models_data %>%
        # Complete missing dates
        tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days")) %>%
        # Creacion de dataset por estacion
        dplyr::mutate(year          = lubridate::year(date),
                      month         = lubridate::month(date),
                      day           = lubridate::day(date),
                      doy           = lubridate::yday(date),
                      time          = as.numeric(date)/1000,
                      season        = lubridate::quarter(date, fiscal_start = 12),
                      prcp_occ_prev = lag(prcp_occ),
                      tipo_dia_prev = lag(tipo_dia),
                      prcp_amt_prev = lag(prcp_amt),
                      tmax_prev     = lag(tmax),
                      tmin_prev     = lag(tmin),
                      row_num       = row_number()) %>%
        # Add seasonal climate to station_climate
        dplyr::left_join(seasonal_climate,
                         by = c("station_id", "year", "season")) %>%
        # Add seasonal covariates to station_climate
        dplyr::mutate(ST1 = if_else(season == 1, seasonal_prcp, 0),
                      ST2 = if_else(season == 2, seasonal_prcp, 0),
                      ST3 = if_else(season == 3, seasonal_prcp, 0),
                      ST4 = if_else(season == 4, seasonal_prcp, 0),
                      SN1 = if_else(season == 1, seasonal_tmin, 0),
                      SN2 = if_else(season == 2, seasonal_tmin, 0),
                      SN3 = if_else(season == 3, seasonal_tmin, 0),
                      SN4 = if_else(season == 4, seasonal_tmin, 0),
                      SX1 = if_else(season == 1, seasonal_tmax, 0),
                      SX2 = if_else(season == 2, seasonal_tmax, 0),
                      SX3 = if_else(season == 3, seasonal_tmax, 0),
                      SX4 = if_else(season == 4, seasonal_tmax, 0))
    tiempo.add_covariates <- proc.time() - t

    # Ungroup and collect data from cluster
    models_data <- models_data %>%
        dplyr::ungroup()

    ##############
    # Progress Bar
    pb$tick(5, tokens = list(c = 0, what = "init fit", stn = "--"))

    # Add stations's latitude and longitude to models_data
    t <- proc.time()
    models_data <- models_data %>%
    dplyr::left_join(stations %>% dplyr::select(station_id, latitude, longitude),
                     by = 'station_id')
    tiempo.join_stations <- proc.time() - t


    ##############
    # Progress Bar
    pb$tick(5, tokens = list(c = 0, what = "init fit", stn = "--"))


    ######################################
    ## Control de tiempo de la preparación
    tiempo.prep_data <- proc.time() - t.p



    #########################
    ## INICIAN LOS AJUSTES ##
    #########################


    ################################################
    # Progress Bar (for parallel execution) ----
    if(nworkers > 1)
        pb$tick(0, tokens = list(c = 0, what = "fitting"))


    t.m <- proc.time()
    models <- foreach::foreach(station = unique_stations, .multicombine = T, .packages = c('dplyr'),
                               .options.snow = list(progress = progress_pb), .verbose = verbose) %dopar% {

        station_climate <- models_data %>%
            dplyr::filter(station_id == station)


        ##########################################
        ## Fit model for precipitation occurrence.

        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(0, tokens = list(c = which(unique_stations == station), what = "prcp_occ", stn = station))

        # Remove NAs
        prcp_occ_fit_noNA_cols <- c('prcp_occ', 'prcp_occ_prev', 'ST1', 'ST2', 'ST3', 'ST4',
                                    'doy', 'time', 'row_num')
        prcp_occ_indexes <- station_climate %>% tidyr::drop_na(prcp_occ_fit_noNA_cols) %>% dplyr::pull(row_num)

        # Create formula
        prcp_occ_fm <- prcp_occ ~ s(tipo_dia_prev, bs = 're') +
            s(time, bs = 'gp', k = 100) + s(doy, bs = 'cc', k = 20)

        if (control$use_covariates) {
            prcp_occ_cov <- models_data %>% dplyr::select(dplyr::matches('ST\\d')) %>% names
            prcp_occ_cov_fm_str <- paste("s(", prcp_occ_cov, ",",
                                         "bs = c('tp' ), k = 20)", collapse = " + ")
            prcp_occ_cov_fm     <- stats::as.formula(paste('~ . +', prcp_occ_cov_fm_str))

            prcp_occ_fm     <- stats::update(prcp_occ_fm, prcp_occ_cov_fm )
        }

        t <- proc.time()
        prcp_occ_fit <- mgcv::bam(formula = prcp_occ_fm,
                                  data = station_climate[prcp_occ_indexes,],
                                  family = stats::binomial(probit),
                                  method = 'fREML',
                                  #cluster = cluster,
                                  control = list(nthreads = control$avbl_cores))
        tiempo.prcp_occ_fit <- proc.time() - t

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        station_climate$prcp_occ_residuals <- NA
        station_climate$prcp_occ_residuals[prcp_occ_indexes] <- residuals(prcp_occ_fit, type = 'response')

        # Se agregan datos de control
        prcp_occ_fit[["fitted_station"]] <- station
        prcp_occ_fit[["dates_used_fitting"]] <- station_climate[prcp_occ_indexes, ] %>%
            dplyr::pull(date)
        prcp_occ_fit[["residuals_tibble"]]   <- station_climate[prcp_occ_indexes, ] %>%
            dplyr::select(date, residual = prcp_occ_residuals)
        prcp_occ_fit[["execution_time"]]     <- tiempo.prcp_occ_fit
        # End of prcp_occ fit

        ## Fit model for precipitation occurrence.
        #########################################


        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(30, tokens = list(c = which(unique_stations == station), what = "prcp_occ", stn = station))

        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(0, tokens = list(c = which(unique_stations == station), what = "prcp_amt", stn = station))


        #######################################
        ## Fit model for precipitation amounts.

        # Remove NAs
        prcp_amt_fit_noNA_cols <- c('prcp_amt', 'prcp_occ_prev', 'ST1', "ST2", 'ST3', 'ST4',
                                    'doy', 'time', 'row_num')
        gamma_indexes <- station_climate %>% tidyr::drop_na(prcp_amt_fit_noNA_cols) %>% dplyr::pull(row_num)

        # Create formula
        prcp_amt_fm <- prcp_amt ~ s(prcp_occ_prev, bs = 're')

        if (control$use_covariates) {
            prcp_amt_cov <- models_data %>% dplyr::select(dplyr::matches('seasonal_prcp')) %>% names
            prcp_amt_cov_str     <- paste("s(", prcp_amt_cov, ",",
                                          "bs = c('tp' ), k = 20)")
            prcp_amt_cov_fm     <- stats::as.formula(paste('~ . +', prcp_amt_cov_str))
            # prcp_amt_cov_fm <- ~ . + ST1 + ST2 + ST3 + ST4
            prcp_amt_fm         <- stats::update( prcp_amt_fm, prcp_amt_cov_fm )
        }


        prcp_amt_fit <- foreach::foreach(m = 1:12, .multicombine = T) %dopar% {
            t <- proc.time()
            prcp_amt_fit <- mgcv::bam(formula = prcp_amt_fm,
                                      data = station_climate[gamma_indexes,] %>%
                                          dplyr::filter(as.logical(prcp_occ) & month == m),
                                      family = stats::Gamma(link = log),
                                      method = 'fREML',
                                      #cluster = cluster,
                                      control = list(nthreads = control$avbl_cores))
            tiempo.prcp_amt_fit <- proc.time() - t

            # Se agregan las fechas para poder hacer las validaciones posteriormente
            prcp_amt_fit[["fitted_station"]] <- station
            # Se agregan datos de control
            prcp_amt_fit[["dates_used_fitting"]] <- station_climate[gamma_indexes, ] %>%
                dplyr::filter(as.logical(prcp_occ) & month == m) %>%
                dplyr::pull(date)
            prcp_amt_fit[["execution_time"]]     <- tiempo.prcp_amt_fit
            # End of prcp_amt fit

            return (prcp_amt_fit)
        }

        ## Fit model for precipitation amounts.
        #######################################


        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(20, tokens = list(c = which(unique_stations == station), what = "prcp_amt", stn = station))

        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(0, tokens = list(c = which(unique_stations == station), what = "tmax", stn = station))


        ################################
        ## Fit model for max temperature.

        # Remove NAs
        tmax_fit_noNA_cols <- c('tmax', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev',
                                'seasonal_tmax', 'seasonal_tmin', 'row_num')
        tmax_indexes <- station_climate %>% tidyr::drop_na(tmax_fit_noNA_cols) %>% dplyr::pull(row_num)
        # Fit model for max temperature.

        # Create formula
        tmax_fm <- tmax ~ s(tmax_prev, tmin_prev, k = 50) +
            s(prcp_occ, bs = "re") +
            s(prcp_occ_prev, bs = "re") +
            s(time, bs = 'gp', k = 100) +
            s(doy, bs = "cc", k = 30)

        if (control$use_covariates) {
            tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
            tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
            tmax_cov_fm_str <- paste("s(", tmax_cov, ", ", tmin_cov, ",",
                                     "k = 20)", collapse = " + ")
            tmax_cov_fm     <- stats::as.formula(paste('~ . +', tmax_cov_fm_str))
            tmax_fm         <- stats::update( tmax_fm, tmax_cov_fm )
        }

        t <- proc.time()
        tmax_fit <- mgcv::bam(formula = tmax_fm,
                              data = station_climate[tmax_indexes,],
                              method = 'fREML',
                              #cluster = cluster,
                              control = list(nthreads = control$avbl_cores))
        tiempo.tmax_fit <- proc.time() - t

        station_climate$tmax_residuals <- NA
        station_climate$tmax_residuals[tmax_indexes] <- residuals(tmax_fit, type = 'response')

        # Se agregan datos de control
        tmax_fit[["fitted_station"]] <- station
        tmax_fit[["dates_used_fitting"]] <- station_climate[tmax_indexes, ] %>%
            dplyr::pull(date)
        tmax_fit[["residuals_tibble"]]   <- station_climate[tmax_indexes, ] %>%
            dplyr::select(date, residual = tmax_residuals)
        tmax_fit[["execution_time"]]     <- tiempo.tmax_fit
        # End of tx fit

        ## Fit model for max temperature.
        ################################


        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(20, tokens = list(c = which(unique_stations == station), what = "tmax", stn = station))

        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(0, tokens = list(c = which(unique_stations == station), what = "tmin", stn = station))


        #################################
        ## Fit model for min temperature.

        # Fit model for min temperature.
        tmin_fit_noNA_cols <- c('tmin', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev',
                                'seasonal_tmax', 'seasonal_tmin', 'row_num')
        tmin_indexes <- station_climate %>% tidyr::drop_na(tmin_fit_noNA_cols) %>% dplyr::pull(row_num)

        # Create formula
        tmin_fm <- tmin ~ s(tmax_prev, tmin_prev, k = 50) +
            s(prcp_occ, bs = "re") +
            s(prcp_occ_prev, bs = "re") +
            s(time, bs = 'gp', k = 100) +
            s(doy, bs = "cc", k = 30)

        if (control$use_covariates) {
            tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
            tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
            tmin_cov_fm_str <- paste("s(", tmax_cov, ", ", tmin_cov, ",",
                                     "k = 20)", collapse = " + ")
            tmin_cov_fm     <- stats::as.formula(paste('~ . +', tmin_cov_fm_str))
            tmin_fm         <- stats::update( tmin_fm, tmin_cov_fm )
        }

        t <- proc.time()
        tmin_fit <- mgcv::bam(formula = tmin_fm,
                              data = station_climate[tmin_indexes,],
                              method = 'fREML',
                              #cluster = cluster,
                              control = list(nthreads = control$avbl_cores))
        tiempo.tmin_fit <- proc.time() - t

        station_climate$tmin_residuals <- NA
        station_climate$tmin_residuals[tmin_indexes] <- residuals(tmin_fit, type = 'response')


        # Se agregan datos de control
        tmin_fit[["fitted_station"]] <- station
        tmin_fit[["dates_used_fitting"]] <- station_climate[tmin_indexes, ] %>%
            dplyr::pull(date)
        tmin_fit[["residuals_tibble"]]   <- station_climate[tmin_indexes, ] %>%
            dplyr::select(date, residual = tmin_residuals)
        tmin_fit[["execution_time"]]     <- tiempo.tmin_fit
        # End of tn fit

        ################################################
        # Progress Bar (for non parallel execution) ----
        if(nworkers == 1)
            pb$tick(20, tokens = list(c = which(unique_stations == station), what = "tmin", stn = station))


        # Residuals estimation
        residuos <- station_climate %>%
            dplyr::select(station_id, date, dplyr::ends_with("residuals"), tipo_dia)

        estadisticas.umbrales <- purrr::map(
            .x = 1:12,
            .f = function(mes) {
                umbrales.mes <- station_climate %>%
                    dplyr::filter(lubridate::month(date) == mes)

                umbrales.con.prcp <- dplyr::filter(umbrales.mes, prcp > 0.1) %>%
                    dplyr::group_by(., lubridate::month(date)) %>%
                    dplyr::summarise(., min.range = min(tmax - tmin, na.rm = T),
                                     max.range = max(tmax - tmin, na.rm = T)) %>%
                    dplyr::select(., min.range, max.range)

                umbrales.sin.prcp <- dplyr::filter(umbrales.mes, prcp <= 0.1) %>%
                    dplyr::group_by(., lubridate::month(date)) %>%
                    dplyr::summarise(., min.range = min(tmax - tmin, na.rm = T),
                                     max.range = max(tmax - tmin, na.rm = T))  %>%
                    dplyr::select(., min.range, max.range)

                return (list(con.prcp = umbrales.con.prcp, sin.prcp = umbrales.sin.prcp))
            }
        )

        # Retornar resultados como lista
        return_list <- list(station_id = station,
                            models_data = station_climate,
                            gam_fits  = list(prcp_occ_fit = prcp_occ_fit,
                                             prcp_amt_fit = prcp_amt_fit,
                                             tmax_fit = tmax_fit,
                                             tmin_fit = tmin_fit),
                            residuals = residuos,
                            estadisticos.umbrales = estadisticas.umbrales)

        return(return_list)

    }
    tiempo.models <- proc.time() - t.m


    ##############
    # Progress Bar
    pb$tick(1, tokens = list(c = length(unique_stations), what = "fit finished", stn = "--"))


    ###########################
    ## Preparar datos de salida

    # Save start_climatology to returned model
    model[["start_climatology"]] <- models_data %>%
        dplyr::select(station_id, date, tmax, tmin, prcp_occ) %>%
        dplyr::group_by(station_id, month = lubridate::month(date), day = lubridate::day(date)) %>%
        dplyr::summarise(tmax = mean(tmax, na.rm = T),
                         tmin = mean(tmin, na.rm = T),
                         prcp_occ = mean(prcp_occ, na.rm = T)) %>%
        dplyr::ungroup()

    # Save gam results to returned model
    gam_fits <- lapply(models, '[[', 'gam_fits')
    names(gam_fits) <- lapply(models, '[[', 'station_id')
    model[["gam_fits"]] <- gam_fits

    # Save residuals to returned model
    residuals <- lapply(models, '[[', 'residuals')
    names(residuals) <- lapply(models, '[[', 'station_id')
    model[["models_residuals"]] <- residuals

    # # Save estadisticos.residuos to returned model
    # estadisticos.residuos <- lapply(models, '[[', 'estadisticos.residuos')
    # names(estadisticos.residuos) <- lapply(models, '[[', 'station_id')
    # model[["estadisticos.residuos"]] <- estadisticos.residuos

    # Save estadisticos.umbrales to returned model
    estadisticos.umbrales <- lapply(models, '[[', 'estadisticos.umbrales')
    names(estadisticos.umbrales) <- lapply(models, '[[', 'station_id')
    model[["estadisticos.umbrales"]] <- estadisticos.umbrales

    # Save execution time to returned model
    model[["execution_time"]] <- tiempo.models

    # Set model's class
    class(model) <- "gamwgen"


    #########################
    ## FINALIZAR EJECUCIÓN ##
    #########################


    ## Cerrar progress bar
    pb$terminate()


    ## Remove parallelization conf, if necessary
    if(remove_parallelization_conf) {
        foreach::registerDoSEQ()
        snow::stopCluster(cluster)
    }

    ## Return model
    model
}
