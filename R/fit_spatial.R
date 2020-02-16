
#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
spatial_fit_control <- function(prcp_occurrence_threshold = 0.1,
                                avbl_cores = 2,
                                planar_crs_in_metric_coords = 22185) {

    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                avbl_cores = avbl_cores,
                planar_crs_in_metric_coords = planar_crs_in_metric_coords))
}


#' @title Weather model calibration
#' @description Fits a weather model from historical data.
#' @import dplyr
#' @export
spatial_calibrate <- function(climate, stations, seasonal_covariates = NULL,
                              control = gamwgen:::spatial_fit_control(),
                              verbose = F) {

    ## Español: Se crea el objeto a ser retornado al culminar el ajuste!
    ## English: An object is created to be returned once the fitting process ends!
    model <- list()

    ###############################################################

    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("foreach"))

    ###############################################################

    # Español: Se controlan que los datos recibidos tengan el formato correcto
    # English: It is checked that the input data is in the correct format
    gamwgen:::check.fit.input.data(climate, stations, seasonal_covariates)

    ###############################################################

    # Español: Se crean crean las variables complementarias que serán utilizadas en el ajuste
    # English: omplementary variables are created that will be used in model fitting
    climate <- climate %>%
        dplyr::arrange(station_id, date) %>%
        dplyr::mutate(prcp_occ_threshold = control$prcp_occurrence_threshold) %>%
        dplyr::mutate(prcp_occ = as.integer(prcp >= control$prcp_occurrence_threshold),  # prcp occurrence
                      type_day = factor(ifelse(as.logical(prcp_occ), 'Wet', 'Dry'), levels = c('Wet', 'Dry')),
                      prcp_amt = ifelse(as.logical(prcp_occ), prcp, NA_real_))  # prcp amount/intensity

    # Español: Se controla que las temperaturas máximas no sean inferiores a las mínimas. Si se detectan,
    # los valores son reemplazados por NAs
    # English: It is checked that maximum temperature are not below minimum temperatures. If detected,
    # values are replaced by NAs
    invalid_records <- c(which(climate$tmax < climate$tmin), which(climate$prcp < 0))

    if(length(invalid_records) > 0) {
        str_stns <- paste(unique(climate[invalid_records, 'station_id']), collapse = ", " )
        warning(glue::glue('Setting to NA {length(invalid_records)} records with maximum ',
                           'temperature < minimum temperature and/or negative rainfalls. ',
                           'The afected stations are: {str_stns}.'))
        climate[invalid_records, c('tmax', 'tmin', 'prcp')] <- NA
    }

    ###############################################################

    # Español: Se convierte el objeto con las coordenadas de las estaciones
    # al sistema de proyección definido
    # English: The object with the coordinates of the weather stations is
    # converted to the defined reference system
    stations <- stations %>%
        sf::st_transform(control$planar_crs_in_metric_coords) %>%
        dplyr::mutate(longitude = sf::st_coordinates(geometry)[,'X'],
                      latitude  = sf::st_coordinates(geometry)[,'Y']) %>%
        dplyr::arrange(station_id)

    # Español: Se obtienen los índices de las estaciones en el data frame
    # con datos climaticos
    # English: The weather stations' index in the climate data frame
    # are extracted
    unique_stations <- sort(unique(climate$station_id))

    # Español: Se verifique que las estaciones en el df con datos climaticos
    # sean las mismas que en el df de estaciones
    # English: it is verified that the weather stations in the climate data frame
    # are the same in the stations data frame
    if (!all(unique_stations == stations$station_id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    ###############################################################

    # Español: Si serán utilizadas, comienza la preparación y verificación del
    # objeto con las covariables estacionales.
    # English: If they are to be used, begins the preparation and verification of
    # the object with seasonal covariates
    if (!is.null(seasonal_covariates))
        seasonal_covariates <- seasonal_covariates %>% dplyr::arrange(station_id, year, season)

    # Español: Se verifica que hayan covariables suficientes para cubrir
    # todas las fechas en climate
    # English: It is verified that the covariates are enough to cover
    # all dates in the climate data frame
    if (!is.null(seasonal_covariates)) {
        climate_control <- climate %>% dplyr::distinct(station_id, year = as.integer(lubridate::year(date)),
                                                       season = lubridate::quarter(date, fiscal_start = 12))
        seasnl_cov_ctrl <- seasonal_covariates %>% dplyr::distinct(station_id, year, season)
        if (!dplyr::all_equal(climate_control, seasnl_cov_ctrl))
            stop("Years in climate and years in seasonal_covariates don't match!")
    }

    ###############################################################



    ######################################
    ## AQUI EMPIEZA REALMENTE EL AJUSTE ##
    ######################################


    #################################
    ## Español: Control de tiempo de ejecución
    ## English: Control of the execution time
    t.m <- proc.time()


    ###########################################
    ## Español: Se guardan en el objeto model los datos de entrada
    ## English: Inpuat data is stored in the object model

    model[["control"]] <- control

    model[["stations"]] <- stations

    model[["climate"]] <- climate

    model[["seasonal_covariates"]] <- seasonal_covariates

    model[["crs_used_to_fit"]] <- sf::st_crs(control$planar_crs_in_metric_coords)



    ####################################
    ## Parallelization initialization ##
    ####################################

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



    ###################################
    ## PREPARACIÓN BARRA DE PROGRESO ##
    ###################################


    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = " fitting: :what | previous progress: [:bar] :percent (in :elapsed)",
        total = 101, clear = FALSE, width= 80, show_after = 0)



    ############################
    ## PREPARACIÓN DEL AJUSTE ##
    ############################


    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "init fit"))


    #####################################
    # Español: Control de tiempo de la preparación
    # English: Check the execution time
    t.p <- proc.time()


    # Español: Se agrupan los datos por estación meteorológica
    # English: Group by stations and conf cluster
    models_data <- climate %>%
        dplyr::group_by(station_id)

    # Español: Crear matriz con todas las estaciones
    # English: Creation of a matrix with every station
    t <- proc.time()
    models_data <- models_data %>%
        # Español: Se completan fechas faltantes
        # English: Complete missing dates
        tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days")) %>%
        # Español: Creacion de dataset por estacion
        # English: Creation of a dataset per station
        dplyr::mutate(year          = lubridate::year(date),
                      month         = lubridate::month(date),
                      day           = lubridate::day(date),
                      doy           = lubridate::yday(date),
                      time          = as.numeric(date)/1000,
                      season        = lubridate::quarter(date, fiscal_start = 12),
                      prcp_occ_prev = lag(prcp_occ),
                      type_day_prev = lag(type_day),
                      prcp_amt_prev = lag(prcp_amt),
                      tmax_prev     = lag(tmax),
                      tmin_prev     = lag(tmin))
    # Español: Agregar las covariables si son necesarias
    # English: Add covariables if necesary
    if (!is.null(seasonal_covariates))
        models_data <- models_data %>%
        # Español: Agregar covariables estacionales a station_climate
        # English: Add seasonal climate to station_climate
            dplyr::left_join(seasonal_covariates,
                             by = c("station_id", "year", "season")) %>%
        # Español: Crear covariables estacionales en el formato necesario
        # English: Create seasonal covariates in the needed format
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
    # Español: Calcular el tiempo total para agregar las covariables
    # English: Compute total time to add covariables
    tiempo.add_covariates <- proc.time() - t

    # Español: Desagrupar y calectar los datos del cluster
    # English: Ungroup and collect data from cluster
    models_data <- models_data %>%
        dplyr::ungroup()

    ##############
    # Progress Bar
    pb$tick(5, tokens = list(what = "init fit"))

    # Español: Agregar la latitud y longitud de las estaciones a models_data
    # English: Add stations's latitude and longitude to models_data
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
    ## Español: Control de tiempo de la preparación
    ## English: Control preparation time
    tiempo.prep_data <- proc.time() - t.p



    #########################
    ## INICIAN LOS AJUSTES ##
    #########################


    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "prcp_occ"))


    ##########################################
    ## Español: Ajustar el modelo de ocurrencia de precipitación.
    ## English: Fit model for precipitation occurrence.

    # Español: Remover NAs
    # English: Remove NAs
    prcp_occ_fit_noNA_cols <- c('prcp_occ', 'prcp_occ_prev', 'doy', 'time', 'row_num')
    if (!is.null(seasonal_covariates))
        prcp_occ_fit_noNA_cols <- append(prcp_occ_fit_noNA_cols, c('ST1', 'ST2', 'ST3', 'ST4'))
    prcp_occ_indexes <- models_data %>%
        tidyr::drop_na(tidyselect::all_of(prcp_occ_fit_noNA_cols)) %>%
        dplyr::pull(row_num)

    # Español: Crear formula
    # English: Create formula
    prcp_occ_fm <- prcp_occ ~ s(time, bs = "gp", k = 1000) +
        te(type_day_prev, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations)) +
        te(longitude, latitude, doy, d = c(2, 1), bs = c("tp", "cc"), k = length(unique_stations))

    # Español: Si es necesario, agregar covariables a la formula
    # English: If necesary, add covariates to the formula
    if (!is.null(seasonal_covariates)) {
        prcp_occ_cov <- models_data %>% dplyr::select(dplyr::matches('ST\\d')) %>% names
        prcp_occ_cov_fm_str <- paste("te(", prcp_occ_cov, ", longitude, latitude, d = c(1, 2), ",
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
    # Español: Ajustar el GAM
    # English: Fit the GAM
    prcp_occ_fit <- mgcv::bam(formula = prcp_occ_fm,
                              data = models_data[prcp_occ_indexes,],
                              family = stats::binomial(probit),
                              method = "fREML",
                              cluster = cluster,
                              control = list(nthreads = control$avbl_cores))
    tiempo.prcp_occ_fit <- proc.time() - t

    # Español: Se agregan las residuos para poder ajustar la componente meterológica
    # English: Residues are added to fit the meteorological component
    models_data[prcp_occ_indexes, "prcp_occ_residuals"] <- residuals(prcp_occ_fit, type = 'response')

    # Español: Se agregan datos de control
    # English: Control data is added
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
    ## Español: Ajustar el modelo de montos de precipitación
    ## English: Fit model for precipitation amounts.

    # Español: Remove NAs
    # English: Remove NAs
    prcp_amt_fit_noNA_cols <- c('prcp_amt', 'prcp_occ_prev', 'doy', 'time', 'row_num')
    if (!is.null(seasonal_covariates))
        prcp_amt_fit_noNA_cols <- append(prcp_amt_fit_noNA_cols, c('ST1', 'ST2', 'ST3', 'ST4'))
    gamma_indexes <- models_data %>%
        tidyr::drop_na(tidyselect::all_of(prcp_amt_fit_noNA_cols)) %>%
        dplyr::pull(row_num)

    # Español: Crear formula
    # English: Create formula
    prcp_amt_fm <- prcp_amt ~ te(type_day_prev, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations))

    # Español: Si es necesario, agregar covariables a la formula
    # English: If necesary, add covariates to the formula
    if (!is.null(seasonal_covariates)) {
        prcp_amt_cov_fm     <- ~ . + te(seasonal_prcp, longitude, latitude, d = c(1, 2), bs = c('tp', 'tp'), k = length(unique_stations))
        prcp_amt_fm         <- stats::update( prcp_amt_fm, prcp_amt_cov_fm )
    }

    t.a <- proc.time()
    # Español: Ajustar el GAM
    # English: Fit the GAM
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

            # Español: Se agregan datos de control
            # English: Control data is added
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
    ## Español: Ajustar el modelo de temperatura máxima
    ## English: Fit model for maximum temperature.

    # Español: Remove NAs
    # English: Remove NAs
    tmax_fit_noNA_cols <- c('tmax', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev', 'row_num')
    if (!is.null(seasonal_covariates))
        tmax_fit_noNA_cols <- append(tmax_fit_noNA_cols, c('seasonal_tmax', 'seasonal_tmin'))
    tmax_indexes <- models_data %>%
        tidyr::drop_na(tidyselect::all_of(tmax_fit_noNA_cols)) %>%
        dplyr::pull(row_num)

    # Español: Crear formula
    # English: Create formula
    tmax_fm <- tmax ~ s(time, bs = 'gp', k = 1000) +
        te(tmax_prev, tmin_prev, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        te(type_day, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations)) +
        te(type_day_prev, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations)) +
        te(doy, longitude, latitude, d = c(1, 2), bs = c('cc', 'tp'), k = length(unique_stations))

    # Español: Si es necesario, agregar covariables a la formula
    # English: If necesary, add covariates to the formula
    if (!is.null(seasonal_covariates)) {
        tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
        tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
        tmax_cov_fm_str <- paste("te(", tmax_cov, ", ", tmin_cov, ", longitude, latitude, d = c(2, 2), ",
                                 "k = length(unique_stations))", collapse = " + ")
        tmax_cov_fm     <- stats::as.formula(paste('~ . +', tmax_cov_fm_str))
        # tmax_cov_fm <- ~ . +
        #     te(SX1, SN1, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX2, SN2, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX3, SN3, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX4, SN4, longitude, latitude, d = c(2, 2), k = length(unique_stations))
        tmax_fm         <- stats::update( tmax_fm, tmax_cov_fm )
    }

    # Español: Ajustar el GAM
    # English: Fit the GAM
    t <- proc.time()
    tmax_fit <- mgcv::bam(formula = tmax_fm,
                          data = models_data[tmax_indexes,],
                          method = "fREML",
                          cluster = cluster,
                          control = list(nthreads = control$avbl_cores))
    tiempo.tmax_fit <- proc.time() - t

    # Español: Se agregan las residuos para poder ajustar la componente meterológica
    # English: Residues are added to fit the meteorological component
    models_data[tmax_indexes, "tmax_residuals"] <- residuals(tmax_fit, type = 'response')

    # Español: Se agregan datos de control
    # English: Control data is added
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
    ## Español: Ajustar el modelo de temperatura mínima.
    ## English: Fit model for min temperature.

    # Español: Remove NAs
    # English: Remove NAs
    tmin_fit_noNA_cols <- c('tmin', 'tmax_prev', 'tmin_prev', 'prcp_occ', 'prcp_occ_prev', 'row_num')
    if (!is.null(seasonal_covariates))
        tmin_fit_noNA_cols <- append(tmin_fit_noNA_cols, c('seasonal_tmax', 'seasonal_tmin'))
    tmin_indexes <- models_data %>%
        tidyr::drop_na(tidyselect::all_of(tmin_fit_noNA_cols)) %>%
        dplyr::pull(row_num)

    # Español: Crear formula
    # English: Create formula
    tmin_fm <- tmin ~ s(time, bs = 'gp', k = 1000) +
        te(tmax_prev, tmin_prev, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        te(type_day, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations)) +
        te(type_day_prev, longitude, latitude, d = c(1, 2), bs = c('re', 'tp'), k = length(unique_stations)) +
        te(doy, longitude, latitude, d = c(1, 2), bs = c('cc', 'tp'), k = length(unique_stations))

    # Español: Si es necesario, agregar covariables a la formula
    # English: If necesary, add covariates to the formula
    if (!is.null(seasonal_covariates)) {
        tmax_cov <- models_data %>% dplyr::select(dplyr::matches('SX\\d')) %>% names
        tmin_cov <- models_data %>% dplyr::select(dplyr::matches('SN\\d')) %>% names
        tmin_cov_fm_str <- paste("te(", tmax_cov, ", ", tmin_cov, ", longitude, latitude, d = c(2, 2), ",
                                 "k = length(unique_stations))", collapse = " + ")
        tmin_cov_fm     <- stats::as.formula(paste('~ . +', tmin_cov_fm_str))
        # tmin_cov_fm <- ~ . +
        #     te(SX1, SN1, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX2, SN2, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX3, SN3, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        #     te(SX4, SN4, longitude, latitude, d = c(2, 2), k = length(unique_stations)) +
        tmin_fm         <- stats::update( tmin_fm, tmin_cov_fm )
    }

    # Español: Ajustar el GAM
    # English: Fit the GAM
    t <- proc.time()
    tmin_fit <- mgcv::bam(formula = tmin_fm,
                          data = models_data[tmin_indexes,],
                          method = "fREML",
                          cluster = cluster,
                          control = list(nthreads = control$avbl_cores))
    tiempo.tmin_fit <- proc.time() - t

    # Español: Se agregan las residuos para poder ajustar la componente meterológica
    # English: Residues are added to fit the meteorological component
    models_data[tmin_indexes, "tmin_residuals"] <- residuals(tmin_fit, type = 'response')

    # Español: Se agregan datos de control
    # English: Control data is added
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

    ##############
    # Progress Bar
    pb$tick(0, tokens = list(what = "retry thresholds"))


    ##############
    # Español: Umbral de reintento
    # English: Retry threshold
    statistics_threshold <- purrr::map_dfr(
        .x = 1:12,
        .f = function(mes) {
            threshold.month <- models_data %>%
                dplyr::filter(lubridate::month(date) == mes)

            threshold_wet <- dplyr::filter(umbrales.mes, prcp > control$prcp_occurrence_threshold) %>%
                dplyr::group_by(., station_id, lubridate::month(date)) %>%
                dplyr::summarise(., min.range = min(tmax - tmin, na.rm = T),
                                 max.range = max(tmax - tmin, na.rm = T)) %>%
                dplyr::mutate(., prcp_occ = 1, month = mes) %>%
                dplyr::select(., station_id, month, prcp_occ, min.range, max.range) %>%
                dplyr::ungroup(.)

            threshold_dry <- dplyr::filter(umbrales.mes, prcp <= control$prcp_occurrence_threshold) %>%
                dplyr::group_by(., station_id, lubridate::month(date)) %>%
                dplyr::summarise(., min.range = min(tmax - tmin, na.rm = T),
                                 max.range = max(tmax - tmin, na.rm = T))  %>%
                dplyr::mutate(., prcp_occ = 0, month = mes) %>%
                dplyr::select(., station_id, month, prcp_occ, min.range, max.range) %>%
                dplyr::ungroup(.)


            threshold <- rbind(threshold_wet, threshold_dry)

            return (threshold)
        }
    )


    #################################
    ## Control de tiempo de ejecución
    tiempo.models <- proc.time() - t.m


    ##############
    # Progress Bar
    pb$tick(1, tokens = list(what = "fit finished"))


    ###########################
    ## Preparar datos de salida

    # Español: Se guardan la climatología de cada estación que será
    # usada para inicializar la simulación
    # English: Save weather station's climatology that will be use to
    # initialize the simulation
    model[["start_climatology"]] <- models_data %>%
        dplyr::select(station_id, date, tmax, tmin, prcp_occ) %>%
        dplyr::group_by(station_id, month = lubridate::month(date), day = lubridate::day(date)) %>%
        dplyr::summarise(tmax = mean(tmax, na.rm = T),
                         tmin = mean(tmin, na.rm = T),
                         prcp_occ = mean(prcp_occ, na.rm = T)) %>%
        dplyr::ungroup()

    # Español: Se guardan resultados del GAM en el objeto model
    # English: Save GAM results to returned model
    model[['fitted_models']][["prcp_occ_fit"]] <- prcp_occ_fit
    model[['fitted_models']][["prcp_amt_fit"]] <- prcp_amt_fit
    model[['fitted_models']][["tmax_fit"]]     <- tmax_fit
    model[['fitted_models']][["tmin_fit"]]     <- tmin_fit

    # Español: Se guardan datos usados en el objeto model
    # English: Save models data to returned model
    model[["models_data"]] <- models_data %>%
        dplyr::select(station_id, date, season, prcp, tmax, tmin, prcp_occ, type_day, prcp_amt,
                      prcp_occ_prev, type_day_prev, prcp_amt_prev, tmax_prev, tmin_prev)

    # Español: Se guardan residuos en el objeto model
    # English: Save residuals to returned model
    model[["models_residuals"]] <- models_data %>%
        dplyr::select(station_id, date, dplyr::ends_with("residuals"), type_day)

    # Español: Se guardan los umbrales en el objeto model
    # English: Save estadisticos.umbrales to returned model
    model[["statistics_threshold"]] <- statistics_threshold

    # Español: Se guardan los tiempos de ejecución en el objeto model
    # English: Save execution time to returned model
    model[['exec_times']][["covs_add_time"]] <- tiempo.add_covariates
    model[['exec_times']][["rown_add_time"]] <- tiempo.add_row_numbers
    model[['exec_times']][["join_stn_time"]] <- tiempo.join_stations
    model[['exec_times']][["prep_dat_time"]] <- tiempo.prep_data
    model[['exec_times']][["pocc_fit_time"]] <- tiempo.prcp_occ_fit
    model[['exec_times']][["pamt_fit_time"]] <- tiempo.prcp_amt_fit
    model[['exec_times']][["tmax_fit_time"]] <- tiempo.tmax_fit
    model[['exec_times']][["tmin_fit_time"]] <- tiempo.tmin_fit
    model[['exec_times']][["exec_tot_time"]] <- tiempo.models

    # Español: Se define la clase del objeto
    # English: Set model's class
    class(model) <- "gamwgen"


    #########################
    ## FINALIZAR EJECUCIÓN ##
    #########################


    ## Cerrar progress bar
    pb$terminate()


     ## Remive luster for mgcv
    if (!is.null(cluster))
        parallel::stopCluster(cluster)

    ## Español: Devolver model
    ## English: Return model
    model
}
