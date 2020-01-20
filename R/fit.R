
#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
glmwgen_fit_control <- function(prcp_occurrence_threshold = 0.1,
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
#' @import dplyr
#' @export
calibrate.glmwgen <- function(climate, stations, seasonal_climate = NULL,
                              control = glmwgen:::glmwgen_fit_control(),
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
    glmwgen:::check.fit.input.data(climate, stations, seasonal_climate)

    ###############################################################

    climate <- climate %>%
        dplyr::arrange(station_id, date)

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

    # summarisea_climate es lo que se recibio como seasonal_climate,
    # es decir, seasonal_climate y seasonal_climate son la misma cosa!!
    summarised_climate <- seasonal_climate

    t <- proc.time()
    # Si no se recibió seasonal_climate como parámetro, entonces se la calcula internamente!!
    if (is.null(summarised_climate) | !control$use_external_seasonal_climate)
        summarised_climate <- glmwgen:::summarise_seasonal_climate(climate, control$climate_missing_threshold)
    tiempo.summarised_climate <- proc.time() - t

    # Se verifica que hayan covariables suficientes para cubrir todas las fechas en climate,
    # pero el control solo se hace si se van a utilizar covariables en el ajuste!!
    climate_control <- climate %>% dplyr::distinct(station_id, year = lubridate::year(date),
                                                   season = lubridate::quarter(date, fiscal_start = 12))
    seasnl_cli_ctrl <- summarised_climate %>% dplyr::distinct(station_id, year, season)
    if (!dplyr::all_equal(climate_control, seasnl_cli_ctrl) & control$use_covariates)
        stop("Years in climate and years in seasonal_climate don't match!")

    ###############################################################



    ######################################
    ## AQUI EMPIEZA REALMENTE EL AJUSTE ##
    ######################################


    ###########################################
    ## Parallelization initialization

    ## OBS:
    ## No se usa foreach con %dopar% ni furrr para calcular
    ## los 12 gam en prcp_amt_fit porque el gam ya utiliza
    ## todos los procesadores disponibles. Usando foreach
    ## y %dopar% con la opción cluster del bam, se ebtienen,
    ## al calcular el resultado para prcp_amt_fit, tiempos
    ## 10 veces peores al acutal (30 vs 350)!!
    ## Además, se usa purrr::map en lugar de foreach %do%
    ## porque purrr:map es más rápido!!

    ## Make Cluster
    cluster <- parallel::makeCluster(control$avbl_cores)


    ###########################################
    ## Se guardan en model los datos de entrada

    model[["control"]] <- control

    model[["stations"]] <- stations

    model[["seasonal_data"]] <- summarised_climate

    model[["crs_used_to_fit"]] <- sf::st_crs(control$planar_crs_in_metric_coords)


    #################################
    ## Control de tiempo de ejecución
    t.m <- proc.time()

    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = " fitting: :what | previous progress: [:bar] :percent (in :elapsed)",
        total = 100, clear = FALSE, width= 80, show_after = 0)


    ############################
    ## PREPARACIÓN DEL AJUSTE ##
    ############################


    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "init fit"))


    #####################################
    # Control de tiempo de la preparación
    t.p <- proc.time()


    # Group by stations and conf cluster
    models_data <- climate %>%
        dplyr::mutate(prcp_occ_threshold = control$prcp_occurrence_threshold) %>%
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
                      prcp_occ      = as.integer(prcp >= prcp_occ_threshold),  # prcp occurrence
                      tipo_dia      = factor(ifelse(as.logical(prcp_occ), 'Lluvioso', 'Seco'), levels = c('Lluvioso', 'Seco')),
                      prcp_amt      = ifelse(as.logical(prcp_occ), prcp, NA_real_),  # prcp amount/intensity
                      prcp_occ_prev = lag(prcp_occ),
                      tipo_dia_prev = lag(tipo_dia),
                      prcp_amt_prev = lag(prcp_amt),
                      tmax_prev     = lag(tmax),
                      tmin_prev     = lag(tmin)) %>%
        # Add seasonal covariates to station_climate
        dplyr::left_join(summarised_climate,
                         by = c("station_id", "year", "season")) %>%
        dplyr::mutate(ST1 = if_else(season == 1, sum_prcp, 0),
                      ST2 = if_else(season == 2, sum_prcp, 0),
                      ST3 = if_else(season == 3, sum_prcp, 0),
                      ST4 = if_else(season == 4, sum_prcp, 0),
                      SN1 = if_else(season == 1, mean_tmin, 0),
                      SN2 = if_else(season == 2, mean_tmin, 0),
                      SN3 = if_else(season == 3, mean_tmin, 0),
                      SN4 = if_else(season == 4, mean_tmin, 0),
                      SX1 = if_else(season == 1, mean_tmax, 0),
                      SX2 = if_else(season == 2, mean_tmax, 0),
                      SX3 = if_else(season == 3, mean_tmax, 0),
                      SX4 = if_else(season == 4, mean_tmax, 0))
    tiempo.add_covariates <- proc.time() - t

    # Ungroup and collect data from cluster
    models_data <- models_data %>%
        dplyr::ungroup()

    ##############
    # Progress Bar
    pb$tick(5, tokens = list(what = "init fit"))

    # Add stations's latitude and longitude to models_data
    t <- proc.time()
    models_data <- models_data %>%
    dplyr::left_join(stations %>% dplyr::select(station_id, latitude, longitude),
                     by = 'station_id')
    tiempo.join_stations <- proc.time() - t

    # Se agrega row_num para diferenciar unívocamente cada registro
    t <- proc.time()
    models_data <- models_data %>%
        dplyr::mutate(row_num = row_number())
    tiempo.add_row_numbers <- proc.time() - t


    ##############
    # Progress Bar
    pb$tick(5, tokens = list(what = "init fit"))


    ######################################
    ## Control de tiempo de la preparación
    tiempo.prep_data <- proc.time() - t.p



    #########################
    ## INICIAN LOS AJUSTES ##
    #########################


    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "prcp_occ"))


    ##########################################
    ## Fit model for precipitation occurrence.

    # Remove NAs
    prcp_occ_fit_noNA_cols <- c('prcp_occ', 'prcp_occ_prev', 'ST1', 'ST2', 'ST3', 'ST4',
                                'doy', 'time', 'row_num')
    prcp_occ_indexes <- models_data %>% tidyr::drop_na(prcp_occ_fit_noNA_cols) %>% dplyr::pull(row_num)

    # Create formula
    prcp_occ_fm <- prcp_occ ~ s(tipo_dia_prev, bs = 're') +
        s(latitude, longitude, bs = 'tp', k = length(unique_stations)) +
        s(time, bs = 'gp', k = 20) + s(doy, bs = 'cc', k = 20) +
        te(latitude, longitude, doy, d = c(2, 1), bs = c('tp', 'cc'), k = length(unique_stations))

    if (control$use_covariates) {
        prcp_occ_cov <- models_data %>% dplyr::select(dplyr::matches('ST\\d')) %>% names
        prcp_occ_cov_fm_str <- paste("te(", prcp_occ_cov, ", latitude, longitude, d = c(1, 2), ",
                                     "bs = c('tp' , 'tp'), k = length(unique_stations))", collapse = " + ")
        prcp_occ_cov_fm     <- stats::as.formula(paste('~ . +', prcp_occ_cov_fm_str))
        # prcp_occ_cov_fm <- ~ . +
        #     te(ST1, longitude, latitude, d =c(1, 2), bs = c(“tp” , “tp”), k = length(unique_stations))
        #     te(ST2, longitude, latitude, d =c(1, 2), bs = c(“tp” , “tp”), k = length(unique_stations))
        #     te(ST3, longitude, latitude, d =c(1, 2), bs = c(“tp” , “tp”), k = length(unique_stations))
        #     te(ST4, longitude, latitude, d =c(1, 2), bs = c(“tp” , “tp”), k = length(unique_stations))
        prcp_occ_fm     <- stats::update( prcp_occ_fm, prcp_occ_cov_fm )
    }

    t <- proc.time()
    prcp_occ_fit <- mgcv::bam(formula = prcp_occ_fm,
                              data = models_data[prcp_occ_indexes,],
                              family = stats::binomial(probit),
                              method = "fREML",
                              cluster = cluster,
                              control = list(nthreads = control$avbl_cores))
    tiempo.prcp_occ_fit <- proc.time() - t

    models_data[prcp_occ_indexes, "prcp_occ_residuals"] <- residuals(prcp_occ_fit, type = 'response')

    # Se agregan datos de control
    prcp_occ_fit[["dates_used_fitting"]] <- models_data[prcp_occ_indexes, ] %>%
        dplyr::pull(date)
    prcp_occ_fit[["residuals_tibble"]]   <- models_data[prcp_occ_indexes, ] %>%
        dplyr::select(date, residual = prcp_occ_residuals)
    prcp_occ_fit[["execution_time"]]     <- tiempo.prcp_occ_fit
    # End of prcp_occ fit

    ## Fit model for precipitation occurrence.
    #########################################


    ##############
    # Progress Bar
    pb$tick(30, tokens = list(what = "prcp_occ"))

    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "prcp_amt"))


    #######################################
    ## Fit model for precipitation amounts.

    # Remove NAs
    prcp_amt_fit_noNA_cols <- c('prcp_amt', 'prcp_occ_prev', 'ST1', "ST2", 'ST3', 'ST4',
                                'doy', 'time', 'row_num')
    gamma_indexes <- models_data %>% tidyr::drop_na(prcp_amt_fit_noNA_cols) %>% dplyr::pull(row_num)

    # Create formula
    prcp_amt_fm <- prcp_amt ~ prcp_occ_prev + s(latitude, longitude, k = length(unique_stations))

    if (control$use_covariates) {
        prcp_amt_cov <- models_data %>% dplyr::select(dplyr::matches('ST\\d')) %>% names
        prcp_amt_cov_fm_str <- paste(prcp_amt_cov, collapse = " + ")
        prcp_amt_cov_fm     <- stats::as.formula(paste('~ . +', prcp_amt_cov_fm_str))
        # prcp_amt_cov_fm <- ~ . + ST1 + ST2 + ST3 + ST4
        prcp_amt_fm         <- stats::update( prcp_amt_fm, prcp_amt_cov_fm )
    }

    t.a <- proc.time()
    prcp_amt_fit <- purrr::map(
        .x = 1:12,
        .f = function(m) {
            t <- proc.time()
            prcp_amt_fit_partial <- mgcv::gam(formula = prcp_amt_fm,
                                              data = models_data[gamma_indexes,] %>%
                                                  dplyr::filter(as.logical(prcp_occ) & month == m),
                                              family = stats::Gamma(link = log),
                                              method = 'REML',
                                              cluster = cluster,
                                              control = list(nthreads = control$avbl_cores))
            tiempo.prcp_amt_fit_partial <- proc.time() - t

            # Se agregan datos de control
            prcp_amt_fit_partial[["dates_used_fitting"]] <- models_data[gamma_indexes, ] %>%
                dplyr::filter(as.logical(prcp_occ) & month == m) %>%
                dplyr::pull(date)
            prcp_amt_fit_partial[["execution_time"]]     <- tiempo.prcp_amt_fit_partial
            # End of prcp_amt fit

            return (prcp_amt_fit_partial)
        }
    ) %>% purrr::set_names(readr::date_names_lang("en")$mon_ab)
    tiempo.prcp_amt_fit <- proc.time() - t.a

    ## Fit model for precipitation amounts.
    #######################################


    ##############
    # Progress Bar
    pb$tick(20, tokens = list(what = "prcp_amt"))

    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "tmax"))


    ################################
    ## Fit model for max temperature.

    # Remove NAs
    tx_fit_noNA_cols <- c('tmax', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev',
                          'mean_tmax', 'mean_tmin', 'row_num')
    tmax_indexes <- models_data %>% tidyr::drop_na(tx_fit_noNA_cols) %>% dplyr::pull(row_num)

    # Create formula
    tmax_fm <- tmax ~ te(tmax_prev, tmin_prev, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        prcp_occ + prcp_occ_prev + s(time, bs = 'gp', k = 10) +
        te(latitude, longitude, doy, bs = c('tp', 'cc'), d = c(2, 1), k = length(unique_stations))

    if (control$use_covariates) {
        tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
        tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
        tmax_cov_fm_str <- paste("te(", tmax_cov, ", ", tmin_cov, ", latitude, longitude, d = c(2, 2), ",
                                 "k = length(unique_stations))", collapse = " + ")
        tmax_cov_fm     <- stats::as.formula(paste('~ . +', tmax_cov_fm_str))
        # tmax_cov_fm <- ~ . +
        #     te(SX1, SN1, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX2, SN2, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX3, SN3, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX4, SN4, latitude, longitude, d = c(2, 2), k = length(unique_stations))
        tmax_fm         <- stats::update( tmax_fm, tmax_cov_fm )
    }

    t <- proc.time()
    tmax_fit <- mgcv::bam(formula = tmax_fm,
                          data = models_data[tmax_indexes,],
                          method = "fREML",
                          cluster = cluster,
                          control = list(nthreads = control$avbl_cores))
    tiempo.tmax_fit <- proc.time() - t

    models_data[tmax_indexes, "tmax_residuals"] <- residuals(tmax_fit, type = 'response')

    # Se agregan datos de control
    tmax_fit[["dates_used_fitting"]] <- models_data[tmax_indexes, ] %>%
        dplyr::pull(date)
    tmax_fit[["residuals_tibble"]]   <- models_data[tmax_indexes, ] %>%
        dplyr::select(date, residual = tmax_residuals)
    tmax_fit[["execution_time"]]     <- tiempo.tmax_fit
    # End of tx fit

    ## Fit model for max temperature.
    ################################


    ##############
    # Progress Bar
    pb$tick(20, tokens = list(what = "tmax"))

    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "tmin"))


    #################################
    ## Fit model for min temperature.

    # Remove NAs
    tn_fit_noNA_cols <- c('tmin', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev',
                          'mean_tmax', 'mean_tmin', 'row_num')
    tmin_indexes <- models_data %>% tidyr::drop_na(tn_fit_noNA_cols) %>% dplyr::pull(row_num)

    # Create formula
    tmin_fm <- tmin ~ te(tmin_prev, tmax_prev, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        prcp_occ + prcp_occ_prev + s(time, bs = 'gp', k = 100) +
        te(latitude, longitude, doy, bs = c('tp', 'cc'), d = c(2, 1), k = length(unique_stations))

    if (control$use_covariates) {
        tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
        tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
        tmin_cov_fm_str <- paste("te(", tmax_cov, ", ", tmin_cov, ", latitude, longitude, d = c(2, 2), ",
                                 "k = length(unique_stations))", collapse = " + ")
        tmin_cov_fm     <- stats::as.formula(paste('~ . +', tmin_cov_fm_str))
        # tmin_cov_fm <- ~ . +
        #     te(SX1, SN1, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX2, SN2, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX3, SN3, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX4, SN4, latitude, longitude, d = c(2, 2), k = length(unique_stations)) +
        tmin_fm         <- stats::update( tmin_fm, tmin_cov_fm )
    }

    t <- proc.time()
    tmin_fit <- mgcv::bam(formula = tmin_fm,
                          data = models_data[tmin_indexes,],
                          method = "fREML",
                          cluster = cluster,
                          control = list(nthreads = control$avbl_cores))
    tiempo.tmin_fit <- proc.time() - t

    models_data[tmin_indexes, "tmin_residuals"] <- residuals(tmin_fit, type = 'response')

    # Se agregan datos de control
    tmin_fit[["dates_used_fitting"]] <- models_data[tmin_indexes, ] %>%
        dplyr::pull(date)
    tmin_fit[["residuals_tibble"]]   <- models_data[tmin_indexes, ] %>%
        dplyr::select(date, residual = tmin_residuals)
    tmin_fit[["execution_time"]]     <- tiempo.tmin_fit
    # End of tn fit

    ## Fit model for min temperature.
    #################################


    ##############
    # Progress Bar
    pb$tick(20, tokens = list(what = "tmin"))


    #################################
    ## Control de tiempo de ejecución
    tiempo.models <- proc.time() - t.m


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
    model[['fitted_models']][["prcp_occ_fit"]] <- prcp_occ_fit
    model[['fitted_models']][["prcp_amt_fit"]] <- prcp_amt_fit
    model[['fitted_models']][["tmax_fit"]]     <- tmax_fit
    model[['fitted_models']][["tmin_fit"]]     <- tmin_fit

    # Save models data to returned model
    model[["models_data"]]      <- models_data %>%
        dplyr::select(station_id, date, season, prcp, tmax, tmin, prcp_occ, tipo_dia, prcp_amt,
                      prcp_occ_prev, tipo_dia_prev, prcp_amt_prev, tmax_prev, tmin_prev)

    # Save residuals to returned model
    model[["models_residuals"]] <- models_data %>%
        dplyr::select(station_id, date, dplyr::ends_with("residuals"), tipo_dia)

    # Save execution time to returned model
    model[['exec_times']][["summ_cli_time"]] <- tiempo.summarised_climate
    model[['exec_times']][["covs_add_time"]] <- tiempo.add_covariates
    model[['exec_times']][["rown_add_time"]] <- tiempo.add_row_numbers
    model[['exec_times']][["join_stn_time"]] <- tiempo.join_stations
    model[['exec_times']][["prep_dat_time"]] <- tiempo.prep_data
    model[['exec_times']][["pocc_fit_time"]] <- tiempo.prcp_occ_fit
    model[['exec_times']][["pamt_fit_time"]] <- tiempo.prcp_amt_fit
    model[['exec_times']][["tmax_fit_time"]] <- tiempo.tmax_fit
    model[['exec_times']][["tmin_fit_time"]] <- tiempo.tmin_fit
    model[['exec_times']][["exec_tot_time"]] <- tiempo.models

    # Set model's class
    class(model) <- "glmwgen"


    #########################
    ## FINALIZAR EJECUCIÓN ##
    #########################


    ## Cerrar progress bar
    pb$terminate()


     ## Remive luster for mgcv
    if (!is.null(cluster))
        parallel::stopCluster(cluster)

    ## Return model
    model
}
