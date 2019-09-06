
#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
glmwgen_fit_control <- function(prcp_occurrence_threshold = 0.1,
                                use_external_seasonal_climate = T,
                                climate_missing_threshold = 0.2,
                                gam_control_nthreads = 2,
                                trim_gam_fit_models = F) {

    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                use_external_seasonal_climate = use_external_seasonal_climate,
                climate_missing_threshold = climate_missing_threshold,
                gam_control_nthreads = gam_control_nthreads,
                trim_gam_fit_models = trim_gam_fit_models))
}


#' @title Weather model calibration
#' @description Fits a weather model from historical data.
#' @import foreach
#' @import dplyr
#' @export
calibrate.glmwgen <- function(climate, stations, seasonal.climate = NULL,
                              control = glmwgen:::glmwgen_fit_control(),
                              verbose = T) {

    # Se crea el objeto a ser retornado al culminar el ajuste!
    model <- list()

    ###############################################################

    if (control$use_external_seasonal_climate & is.null(seasonal.climate))
        stop("If use_external_seasonal_climate is True, seasonal.climate can't be NULL!!")

    if (!control$use_external_seasonal_climate & !is.null(seasonal.climate)) {
        warning('Entered seasonal.climate will be descarted because use_external_seasonal_climate is set as False!')
        seasonal.climate <- NULL
    }

    glmwgen:::check.input.data(climate, stations, seasonal.climate)

    ###############################################################

    climate <- climate %>%
        dplyr::rename(station = station_id) %>%
        dplyr::arrange(station, date)

    invalid_records <- c(which(climate$tx < climate$tn), which(climate$prcp < 0))

    if(length(invalid_records) > 0) {
        warning(sprintf('Removing %d records with maximum temperature < minimum temperature or negative rainfalls.', length(invalid_records)))
        climate[invalid_records, c('tx', 'tn', 'prcp')] <- NA
    }

    ###############################################################

    stations <- stations %>%
        dplyr::rename(id = station_id) %>%
        sf::as_Spatial()
    stations <- stations[order(stations$id), ]
    rownames(stations@data) <- stations$id

    # Se verifica si se va a ajustar una o más estaciones.
    fit_multiple_stations <- fit_only_one_station <- F
    if (nrow(stations) == 1) fit_only_one_station <- T
    if (nrow(stations) > 1) fit_multiple_stations <- T

    # Se establecen los nombres de las filas y columnas.
    gral.rownames <- as.character(stations$id)
    gral.colnames <- as.character(stations$id)

    # Se obtienen los estaciones en el df con datos climaticos
    unique_stations <- sort(unique(climate$station))

    # Se verifique que las estaciones en el df con datos climaticos sean las mismas que en el df de estaciones
    if (!all(unique_stations == stations$id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    ###############################################################

    if (control$use_external_seasonal_climate) {

        seasonal.climate <- seasonal.climate %>%
            dplyr::rename(station = station_id)
        seasonal.climate <- seasonal.climate %>%
            dplyr::arrange(station, year, season)

        years_in_climate <- dplyr::distinct(climate, lubridate::year(date)) %>% dplyr::pull()
        years_in_seasonal.climate <- dplyr::distinct(seasonal.climate, year) %>% dplyr::pull()
        if (!dplyr::all_equal(years_in_climate, years_in_seasonal.climate))
            stop("Years in climate and years in seasonal.climate don't match!")
    }

    summarised_climate <- seasonal.climate
    if (is.null(summarised_climate) | !control$use_external_seasonal_climate)
        summarised_climate <- glmwgen:::summarise_seasonal_climate(climate, control$climate_missing_threshold)

    ###############################################################

    ######################################
    ## AQUI EMPIEZA REALMENTE EL AJUSTE ##
    ######################################

    ## TODO: check control variables.
    model[["control"]] <- control

    model[["seasonal"]] <- as.data.frame(summarised_climate)
    model[['stations_proj4string']] <- stations@proj4string
    model[["stations"]] <- stations

    # En caso que se haga el ajuste para múltiples estaciones,
    # es necesario definir una matriz de distancia entre ellas
    if (fit_multiple_stations) {
        dist.mat <- sp::spDists(stations) # no se debe usar al argumento longlat porque fit hace diferentes proyecciones!
        diag(dist.mat) <- 0 # diagonal element is not explicitly equal to zero, so define as such
        colnames(dist.mat) <- gral.colnames
        rownames(dist.mat) <- gral.rownames
        model[["distance_matrix"]] <- dist.mat
    }

    model[["start_climatology"]] <- climate %>%
        group_by(station, month = lubridate::month(date), day = lubridate::day(date)) %>%
        summarise(tx = mean(tx, na.rm = T),
                  tn = mean(tn, na.rm = T),
                  prcp = mean(prcp, na.rm = T))

    # Register a sequential backend if the user didn't register a parallel
    # in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        foreach::registerDoSEQ()
    }

    t.m <- proc.time()
    models <- foreach::foreach(station_id = unique_stations, .multicombine = T, .packages = c('dplyr')) %dopar% {

        station_climate <- climate %>%
            dplyr::filter(station == station_id) %>%
            tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days")) %>%
            dplyr::mutate(year          = lubridate::year(date),
                          month         = lubridate::month(date),
                          day           = lubridate::day(date),
                          doy           = lubridate::yday(date),
                          season        = lubridate::quarter(date, fiscal_start = 12),
                          Rt            = seq(from = -1, to = 1, length.out = length(date)),
                          prcp_occ      = as.integer(prcp >= control$prcp_occurrence_threshold),  # prcp occurrence
                          prcp_amt      = ifelse(as.logical(prcp_occ), prcp, NA_real_),  # prcp amount/intensity
                          prcp_occ_prev = lag(prcp_occ),
                          tx_prev       = lag(tx),
                          tn_prev       = lag(tn),
                          row_num       = row_number()) %>%
            dplyr::left_join(summarised_climate, by = c("station", "year", "season"))


        # Estimate seasonal covariates
        station_summarised_climate  <- summarised_climate %>% dplyr::filter(station == station_id)
        station_seasonal_covariates <- glmwgen:::create_seasonal_covariates(station_summarised_climate)

        # Add seasonal covariates to station_climate
        station_climate <- station_climate %>%
            dplyr::mutate(ST1 = station_seasonal_covariates$prcp[[1]],
                          ST2 = station_seasonal_covariates$prcp[[2]],
                          ST3 = station_seasonal_covariates$prcp[[3]],
                          ST4 = station_seasonal_covariates$prcp[[4]],
                          SMX1 = station_seasonal_covariates$tx[[1]],
                          SMX2 = station_seasonal_covariates$tx[[2]],
                          SMX3 = station_seasonal_covariates$tx[[3]],
                          SMX4 = station_seasonal_covariates$tx[[4]],
                          SMN1 = station_seasonal_covariates$tn[[1]],
                          SMN2 = station_seasonal_covariates$tn[[2]],
                          SMN3 = station_seasonal_covariates$tn[[3]],
                          SMN4 = station_seasonal_covariates$tn[[4]])


        # Fit model for precipitation occurrence.
        prcp_occ_fit_noNA_cols <- c('prcp_occ', 'prcp_occ_prev', 'ST1', 'ST2', 'ST3', 'ST4')
        # station_climate %>% tidyr::drop_na(prcp_occ_fit_noNA_cols) %>% dplyr::summarise_all(funs(sum(is.na(.))))
        probit_indexes <- station_climate %>% tidyr::drop_na(prcp_occ_fit_noNA_cols) %>% dplyr::pull(row_num)

        t <- proc.time()
        prcp_occ_fit <- mgcv::gam(formula = prcp_occ ~ prcp_occ_prev + ST1 + ST2 + ST3 + ST4 + s(doy, bs = "cc"),
                                  select  = TRUE,
                                  data    = station_climate[probit_indexes,],
                                  family  = stats::binomial(probit),
                                  method  = "REML")
        tiempo.prcp_occ_fit <- proc.time() - t

        if (control$trim_gam_fit_models) prcp_occ_fit <- glmwgen:::stripGlmLR(prcp_occ_fit)

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        prcp_occ_fit[["fitted_station"]] <- station_id
        prcp_occ_fit[["fitted_dates"]]   <- station_climate[probit_indexes, ] %>% dplyr::pull(date)
        prcp_occ_fit[["execution_time"]] <- tiempo.prcp_occ_fit
        # End of prcp_occ fit


        # Fit model for precipitation amounts.
        prcp_amt_fit_noNA_cols <- c('prcp_amt', 'prcp_occ_prev', 'ST1', "ST2", 'ST3', 'ST4')
        # station_climate %>% tidyr::drop_na(prcp_amt_fit_noNA_cols) %>% dplyr::summarise_all(funs(sum(is.na(.))))
        gamma_indexes <- station_climate %>% tidyr::drop_na(prcp_amt_fit_noNA_cols) %>% dplyr::pull(row_num)

        prcp_amt_fit <- foreach::foreach(m = 1:12, .multicombine = T) %dopar% {
            t <- proc.time()
            prcp_amt_fit <- mgcv::gam(formula = prcp_amt ~ prcp_occ_prev + ST1 + ST2 + ST3 + ST4,
                                      data    = station_climate[gamma_indexes,] %>% dplyr::filter(month == m),
                                      family  = stats::Gamma(link = log),
                                      method  = 'REML')
            tiempo.prcp_amt_fit <- proc.time() - t

            if (control$trim_gam_fit_models) prcp_amt_fit <- glmwgen:::stripGlmLR(prcp_amt_fit)

            # Se agregan las fechas para poder hacer las validaciones posteriormente
            prcp_amt_fit[["fitted_station"]] <- station_id
            prcp_amt_fit[["fitted_dates"]]   <- station_climate[gamma_indexes, ] %>% dplyr::filter(month == m) %>% dplyr::pull(date)
            prcp_amt_fit[["execution_time"]] <- tiempo.prcp_amt_fit
            # End of prcp_amt fit

            return (prcp_amt_fit)
        }


        # Fit model for max temperature.
        tx_fit_noNA_cols <- c('tx', 'tx_prev', 'tn_prev', 'prcp_occ', 'prcp_occ_prev', 'mean_tx', 'mean_tn')
        # station_climate %>% tidyr::drop_na(tx_fit_noNA_cols) %>% dplyr::summarise_all(funs(sum(is.na(.))))
        tx_indexes <- station_climate %>% tidyr::drop_na(tx_fit_noNA_cols) %>% dplyr::pull(row_num)

        t <- proc.time()
        tx_fit <- mgcv::gam(tx ~ s(tx_prev, tn_prev) + prcp_occ + prcp_occ_prev + s(mean_tx, mean_tn, bs = 'ad') + te(year, doy, bs = 'cc'),
                            select      = TRUE,
                            correlation = nlme::corARMA(form = ~ 1|year, p = 1, q = 2),
                            data        = station_climate[tx_indexes,],
                            method      = "REML",
                            control     = list(nthreads = control$gam_control_nthreads))
        tiempo.tx_fit <- proc.time() - t

        station_climate[tx_indexes, "tx_residuals"] <- tx_fit$residuals

        if (control$trim_gam_fit_models) tx_fit <- glmwgen:::stripGlmLR(tx_fit)

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        tx_fit[["fitted_station"]] <- station_id
        tx_fit[["fitted_dates"]]   <- station_climate[tx_indexes, ] %>% dplyr::pull(date)
        tx_fit[["execution_time"]] <- tiempo.tx_fit
        tx_fit[["residuals.data.frame"]] <- station_climate[tx_indexes, ] %>% dplyr::select(date, residual = tx_residuals)
        # End of tx fit


        # Fit model for min temperature.
        tn_fit_noNA_cols <- c('tn', 'tx_prev', 'tn_prev', 'prcp_occ', 'prcp_occ_prev', 'mean_tx', 'mean_tn')
        # station_climate %>% tidyr::drop_na(tn_fit_noNA_cols) %>% dplyr::summarise_all(funs(sum(is.na(.))))
        tn_indexes <- station_climate %>% tidyr::drop_na(tn_fit_noNA_cols) %>% dplyr::pull(row_num)

        t <- proc.time()
        tn_fit <- mgcv::gam(tn ~ s(tn_prev, tx_prev) + prcp_occ + prcp_occ_prev + s(mean_tn, mean_tx, bs = 'ad') + te(year, doy, bs = 'cc'),
                            select      = TRUE,
                            correlation = corARMA(form = ~ 1|year, p = 1, q = 2),
                            data        = station_climate[tn_indexes,],
                            method      = "REML",
                            control     = list(nthreads = control$gam_control_nthreads))
        tiempo.tn_fit <- proc.time() - t

        station_climate[tn_indexes, "tn_residuals"] <- tn_fit$residuals

        if (control$trim_gam_fit_models) tn_fit <- glmwgen:::stripGlmLR(tn_fit)

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        tn_fit[["fitted_station"]] <- station_id
        tn_fit[["fitted_dates"]]   <- station_climate[tn_indexes, ] %>% dplyr::pull(date)
        tn_fit[["execution_time"]] <- tiempo.tn_fit
        tn_fit[["residuals.data.frame"]] <- station_climate[tn_indexes, ] %>% dplyr::select(date, residual = tn_residuals)
        # End of tn fit


        # Residuals estimation
        residuals_df <- station_climate %>%
            dplyr::select(station, date, month, prcp, prcp_occ, prcp_amt, tx, tn, tx_residuals, tn_residuals) %>%
            dplyr::group_by(month) %>%
            dplyr::mutate(sd.tx_residuals = stats::sd(tx_residuals, na.rm = T),
                          sd.tn_residuals = stats::sd(tn_residuals, na.rm = T)) %>%
            dplyr::group_by(prcp_occ, add = T) %>%
            dplyr::mutate(mean.tx_residuals = mean(tx_residuals, na.rm = T),
                          mean.tn_residuals = mean(tn_residuals, na.rm = T),
                          cov.residuals     = stats::cov(tx_residuals, tn_residuals, use = "pairwise.complete.obs"),
                          var.tx_residuals  = stats::var(tx_residuals, na.rm = T),
                          var.tn_residuals  = stats::var(tn_residuals, na.rm = T)) %>%
            dplyr::mutate(min.range = ifelse(all(is.na(tx-tn)), NA, min(tx - tn, na.rm = T)),
                          max.range = ifelse(all(is.na(tx-tn)), NA, max(tx - tn, na.rm = T))) %>%
            dplyr::ungroup()

        # Crear lista con matrices de covariancia para la creación de ruidos correlacionados
        # mensuales. Se crean 12 matrices para días lluviosos y secos.
        estadisticos.residuos <- foreach::foreach(mes = 1:12, .multicombine = T) %dopar% {

            residuos.mes <- residuals_df %>%
                dplyr::filter(month == mes)

            estadisticos.con.prcp <- residuos.mes %>%
                dplyr::filter(as.logical(prcp_occ)) %>%
                dplyr::distinct(mean.tx_residuals, mean.tn_residuals,
                                cov.residuals,
                                var.tx_residuals, var.tn_residuals) %>%
                list(., est.media = c(.$mean.tx_residuals, .$mean.tn_residuals),
                     est.mat.cov  = matrix(c(.$var.tx_residuals, .$cov.residuals, .$cov.residuals, .$var.tn_residuals), 2, 2)) %>%
                rlang::set_names(c("estadisticos", "media", "matriz.covarianza"))

            estadisticos.sin.prcp <- residuos.mes %>%
                dplyr::filter(!as.logical(prcp_occ)) %>%
                dplyr::distinct(mean.tx_residuals, mean.tn_residuals,
                                cov.residuals,
                                var.tx_residuals, var.tn_residuals) %>%
                list(., est.media = c(.$mean.tx_residuals, .$mean.tn_residuals),
                     est.mat.cov  = matrix(c(.$var.tx_residuals, .$cov.residuals, .$cov.residuals, .$var.tn_residuals), 2, 2)) %>%
                rlang::set_names(c("estadisticos", "media", "matriz.covarianza"))

            return (list(con.prcp = estadisticos.sin.prcp, sin.prcp = estadisticos.sin.prcp))
        }

        # Crear lista con umbrales máximos y mínimos para cada mes.
        estadisticos.umbrales <- foreach::foreach(mes = 1:12, .multicombine = T) %dopar% {

            umbrales.mes <- residuals_df %>%
                dplyr::filter(month == mes)

            umbrales.con.prcp <- umbrales.mes %>%
                dplyr::filter(as.logical(prcp_occ)) %>%
                dplyr::distinct(min.range, max.range)

            umbrales.sin.prcp <- umbrales.mes %>%
                dplyr::filter(!as.logical(prcp_occ)) %>%
                dplyr::distinct(min.range, max.range)

            return (list(con.prcp = umbrales.con.prcp, sin.prcp = umbrales.sin.prcp))
        }


        # Retornar resultados como lista
        return_list <- list(station = station_id,
                            gam_fits  = list(prcp_occ_fit = prcp_occ_fit,
                                             prcp_amt_fit = prcp_amt_fit,
                                             tx_fit = tx_fit,
                                             tn_fit = tn_fit),
                            residuals = dplyr::select(residuals_df, station, date,
                                                      tx_residuals, tn_residuals,
                                                      sd.tx_residuals, sd.tn_residuals),
                            estadisticos.residuos = estadisticos.residuos,
                            estadisticos.umbrales = estadisticos.umbrales)

        return(return_list)
    }
    tiempo.models <- proc.time() - t.m

    # Save gam results to returned model
    gam_fits <- lapply(models, '[[', 'gam_fits')
    names(gam_fits) <- lapply(models, '[[', 'station')
    model[["gam_fits"]] <- gam_fits

    # Save residuals to returned model
    residuals <- lapply(models, '[[', 'residuals')
    names(residuals) <- lapply(models, '[[', 'station')
    model[["residuals"]] <- residuals

    # Save estadisticos.residuos to returned model
    estadisticos.residuos <- lapply(models, '[[', 'estadisticos.residuos')
    names(estadisticos.residuos) <- lapply(models, '[[', 'station')
    model[["estadisticos.residuos"]] <- estadisticos.residuos

    # Save estadisticos.umbrales to returned model
    estadisticos.umbrales <- lapply(models, '[[', 'estadisticos.umbrales')
    names(estadisticos.umbrales) <- lapply(models, '[[', 'station')
    model[["estadisticos.umbrales"]] <- estadisticos.umbrales

    # Save execution time to returned model
    model[["execution_time"]] <- tiempo.models

    class(model) <- "glmwgen"

    model
}
