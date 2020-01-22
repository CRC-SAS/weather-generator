

fit_control <- function(prcp_occurrence_threshold = 0.1,
                                use_seasonal_covariates_precipitation = F, # el modelo más simple es sin covariables
                                use_seasonal_covariates_temperature = F,
                                use_linear_term = F,
                                prcp_lags_to_use = 1,
                                use_amounts = F,
                                tn_lags_to_use = 1,
                                tx_lags_to_use = 1,
                                use_six_months_precipitation = F,
                                use_six_months_temperature = F,
                                save_residuals = F,
                                cov_mode = 'complete',
                                save_lm_fits = F,
                                use_stepwise = F,
                                use_robust_methods = F,
                                use_external_seasonal_climate = T,
                                climate_missing_threshold = 0.2) {

    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                use_seasonal_covariates_precipitation = use_seasonal_covariates_precipitation,
                use_seasonal_covariates_temperature = use_seasonal_covariates_temperature,
                use_linear_term = use_linear_term,
                prcp_lags_to_use = prcp_lags_to_use,
                use_amounts = use_amounts,
                tn_lags_to_use = tn_lags_to_use,
                tx_lags_to_use = tx_lags_to_use,
                use_six_months_precipitation = use_six_months_precipitation,
                use_six_months_temperature = use_six_months_temperature,
                save_residuals = save_residuals,
                cov_mode = cov_mode,
                save_lm_fits = save_lm_fits,
                use_stepwise = use_stepwise,
                use_robust_methods = use_robust_methods,
                use_external_seasonal_climate = use_external_seasonal_climate,
                climate_missing_threshold = climate_missing_threshold))
}


old_fit <- function(climate, stations, seasonal.climate = NULL,
                              control = fit_control(),
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

    unique_stations <- sort(unique(climate$station))

    if (!all(unique_stations == stations$id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    ###############################################################

    # En caso que se haga el ajuste para múltiples estaciones,
    # es necesario definir más variables iniciales!
    if (fit_multiple_stations) {
        dist.mat <- sp::spDists(stations) # no se debe usar al argumento longlat porque fit hace diferentes proyecciones!
        colnames(dist.mat) <- gral.colnames
        rownames(dist.mat) <- gral.rownames
        # diagonal element is not explicitly equal to zero, so define as such
        diag(dist.mat) <- 0
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

    ## TODO: check control variables.
    model[["control"]] <- control

    if (fit_multiple_stations) model[["distance_matrix"]] <- dist.mat
    model[["seasonal"]] <- as.data.frame(summarised_climate)
    model[['stations_proj4string']] <- stations@proj4string
    model[["stations"]] <- stations

    model[["start_climatology"]] <- climate %>%
        group_by(station, month = lubridate::month(date), day = lubridate::day(date)) %>%
        summarise(tx = mean(tx, na.rm = T),
                  tn = mean(tn, na.rm = T),
                  prcp = mean(prcp, na.rm = T))

    # Precipitation covariates
    prcp_covariates <- c("ct", "st")

    # Use six months cycle
    if (control$use_six_months_precipitation) {
        prcp_covariates <- c(prcp_covariates, "ct.six", "st.six")
    }

    if (control$prcp_lags_to_use > 1) {
        prcp_prev_autocorrelation <- c()
        for (i in 2:control$prcp_lags_to_use) {
            prcp_prev_autocorrelation <- c(prcp_prev_autocorrelation,
                                           paste0("prcp_occ_prev.minus.", i))
        }
        prcp_covariates <-c(prcp_covariates, prcp_prev_autocorrelation)
    }

    # Temperature covariates
    temps_covariates <- c("tn_prev", "tx_prev", "ct", "st")

    # Use six months cycle
    if (control$use_six_months_temperature) {
        temps_covariates <- c(temps_covariates, "ct.six", "st.six")
    }

    # Use occurence or amounts
    if (control$use_amounts) {
        temps_covariates <- c(temps_covariates, "prcp")
    } else {
        temps_covariates <- c(temps_covariates, "prcp_occ")
    }

    # Min temperature of n days before day i
    if (control$tn_lags_to_use > 1) {
        tn_prev_autocorrelation <- c()
        for (i in 2:control$tn_lags_to_use) {
            tn_prev_autocorrelation <- c(tn_prev_autocorrelation,
                                         paste0("tn_prev.minus.", i))
        }
        temps_covariates <-c(temps_covariates, tn_prev_autocorrelation)
    }

    # Max temperature of n days before day i
    if (control$tx_lags_to_use > 1) {
        tx_prev_autocorrelation <- c()
        for (i in 2:control$tx_lags_to_use) {
            tx_prev_autocorrelation <- c(tx_prev_autocorrelation,
                                         paste0("tx_prev.minus.", i))
        }
        temps_covariates <-c(temps_covariates, tx_prev_autocorrelation)
    }

    # Use linear term
    if(control$use_linear_term) {
        temps_covariates <- c(temps_covariates, "Rt")
    }

    # Use seasonal covariates
    if (control$use_seasonal_covariates_precipitation) {
        prcp_covariates <- c(prcp_covariates, "ST1", "ST2", "ST3", "ST4")
    }

    if (control$use_seasonal_covariates_temperature) {
        temps_covariates <- c(temps_covariates, "SMN1", "SMN2", "SMN3", "SMN4", "SMX1", "SMX2", "SMX3", "SMX4")
    }

    # Register a sequential backend if the user didn't register a parallel
    # in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        foreach::registerDoSEQ()
    }

    full_dates <- data.frame(date = seq.Date(min(climate$date), max(climate$date), by = 'days'))

    models <- foreach::foreach(station_id = unique_stations, .multicombine = T, .packages = c('dplyr')) %dopar% {
        station_climate <- full_dates %>% left_join(climate[climate$station == station_id, ], by = 'date') %>%
            mutate(prcp_occ = (prcp >= control$prcp_occurrence_threshold) + 0,              # prcp occurrence
                   prcp_amt = ifelse(prcp >= control$prcp_occurrence_threshold, prcp, NA),  # prcp amount/intensity
                   prcp_occ_prev = lag(prcp_occ),
                   tx_prev = lag(tx),
                   tn_prev = lag(tn),
                   tmean = (tx + tn)/2,
                   year_fraction = 2 * pi * lubridate::yday(date)/ifelse(lubridate::leap_year(date), 366, 365),
                   ct = cos(year_fraction),
                   st = sin(year_fraction),
                   Rt = seq(from = -1, to = 1, length.out = length(year_fraction)),
                   row_num = row_number())

        station_climate$station <- station_id

        # Add six months cycle functions
        if (control$use_six_months_precipitation | control$use_six_months_temperature) {
            year_fraction = 2 * pi * lubridate::yday(station_climate$date)/ifelse(lubridate::leap_year(station_climate$date), 366/2, 365/2)
            ct.six = cos(year_fraction)
            st.six = sin(year_fraction)

            station_climate <- data.frame(station_climate,
                                          ct.six,
                                          st.six)
        }

        # Add lagged precipitation
        if (control$prcp_lags_to_use > 1) {
            prcp_occ_prev.mas.i <- as.data.frame(matrix(ncol = control$prcp_lags_to_use-1, nrow = nrow(station_climate)))
            for (i in 2:control$prcp_lags_to_use) {
                prcp_prev <- lag(station_climate$prcp_occ, i)
                prcp_occ_prev.mas.i[,i-1] <- prcp_prev
            }

            colnames(prcp_occ_prev.mas.i) <- prcp_prev_autocorrelation
            station_climate <- data.frame(station_climate,
                                          prcp_occ_prev.mas.i)
        }

        # Add lagged min temperature
        if (control$tn_lags_to_use > 1) {
            tn_prev.mas.i <- as.data.frame(matrix(ncol = control$tn_lags_to_use-1, nrow = nrow(station_climate)))
            for (i in 2:control$tn_lags_to_use) {
                tn_prev <- lag(station_climate$tn, i)
                tn_prev.mas.i[,i-1] <- tn_prev
            }

            colnames(tn_prev.mas.i) <- tn_prev_autocorrelation
            station_climate <- data.frame(station_climate,
                                          tn_prev.mas.i)
        }

        # Add lagged max temperature
        if (control$tx_lags_to_use > 1) {
            tx_prev.mas.i <- as.data.frame(matrix(ncol = control$tx_lags_to_use-1, nrow = nrow(station_climate)))
            for (i in 2:control$tx_lags_to_use) {
                tx_prev <- lag(station_climate$tx, i)
                tx_prev.mas.i[,i-1] <- tx_prev
            }

            colnames(tx_prev.mas.i) <- tx_prev_autocorrelation
            station_climate <- data.frame(station_climate,
                                          tx_prev.mas.i)
        }

        # Estimate seasonal covariates
        if (control$use_seasonal_covariates_precipitation | control$use_seasonal_covariates_temperature) {
            station_summarised_climate  <- summarised_climate %>% dplyr::filter(station == station_id)
            station_seasonal_covariates <- glmwgen:::create_seasonal_covariates(station_summarised_climate)

            station_climate <- data.frame(station_climate)

            if (control$use_seasonal_covariates_precipitation) {
                station_climate <- station_climate %>%
                    dplyr::mutate(ST1 = station_seasonal_covariates$prcp[[1]],
                                  ST2 = station_seasonal_covariates$prcp[[2]],
                                  ST3 = station_seasonal_covariates$prcp[[3]],
                                  ST4 = station_seasonal_covariates$prcp[[4]])
            }
            if (control$use_seasonal_covariates_temperature) {
                station_climate <- station_climate %>%
                    dplyr::mutate(SMX1 = station_seasonal_covariates$tx[[1]],
                                  SMX2 = station_seasonal_covariates$tx[[2]],
                                  SMX3 = station_seasonal_covariates$tx[[3]],
                                  SMX4 = station_seasonal_covariates$tx[[4]],
                                  SMN1 = station_seasonal_covariates$tn[[1]],
                                  SMN2 = station_seasonal_covariates$tn[[2]],
                                  SMN3 = station_seasonal_covariates$tn[[3]],
                                  SMN4 = station_seasonal_covariates$tn[[4]])
            }
        }

        # Fit model for precipitation occurrence.
        probit_indexes <- na.omit(station_climate[, c("prcp_occ", "row_num", "prcp_occ_prev", prcp_covariates)])$row_num

        if (control$use_robust_methods) {
            prcp_occ_fit <- robust::glmRob(formula(paste0("prcp_occ", "~", "prcp_occ_prev +", paste0(c(prcp_covariates), collapse = "+"))),
                                           data = station_climate[probit_indexes,],
                                           family = stats::binomial(probit))
        } else {
            prcp_occ_fit <- stats::glm(formula(paste0("prcp_occ", "~", "prcp_occ_prev +", paste0(c(prcp_covariates), collapse = "+"))),
                                       data = station_climate[probit_indexes, ],
                                       family = stats::binomial(probit))
        }
        coefocc <- coefficients(prcp_occ_fit)

        if(control$use_stepwise) {
            prcp_occ_fit <- stats::step(prcp_occ_fit, trace = 0)
            coefocc[names(coefocc)] <- 0
            coefocc[names(coefficients(prcp_occ_fit))] <- coefficients(prcp_occ_fit)
        }

        station_climate[probit_indexes, "probit_residuals"] <- prcp_occ_fit$residuals

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        prcp_occ_fit[["dates"]] <- station_climate[probit_indexes, "date"]
        # End of prcp_occ fit


        # Fit model for precipitation amounts.
        gamma_indexes <- na.omit(station_climate[, c("prcp_amt", "row_num", prcp_covariates)])$row_num
        # save model in list element 'i'
        if (control$use_robust_methods) {
            prcp_amt_fit <- robustbase::glmrob(formula(paste0("prcp_amt", "~", paste0(prcp_covariates, collapse = "+"))),
                                               data = station_climate[gamma_indexes, ],
                                               family = stats::Gamma(link = 'log'),
                                               method="Mqle")
            alphamt  <- gamma.shape.glm.rob(prcp_amt_fit)$alpha
        } else {
            prcp_amt_fit <- stats::glm(formula(paste0("prcp_amt", "~", paste0(prcp_covariates, collapse = "+"))),
                                       data = station_climate[gamma_indexes, ],
                                       family = stats::Gamma(link = log))
            alphamt  <- MASS::gamma.shape(prcp_amt_fit)$alpha
        }
        coefamt <- prcp_amt_fit$coefficients

        # if(control$use_stepwise) {
        #     prcp_amt_fit <- stats::step(prcp_amt_fit)
        # }

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        prcp_amt_fit[["dates"]] <- station_climate[gamma_indexes, "date"]
        # End of prcp_amt fit


        # Fit model for max temperature.
        tx_indexes <- na.omit(station_climate[, c("tx", "row_num", temps_covariates)])$row_num
        # como se agrega tn al ajuste de tx, entonces, también se debe excluir los NAs en tn
        tx_indexes <- na.omit(station_climate[tx_indexes, c("tn", "row_num", temps_covariates)])$row_num

        if (control$use_robust_methods) {
            tx_fit <- robustbase::lmrob(formula(paste0("tx", "~", "tn", "+", paste0(temps_covariates, collapse = "+"))),
                                        data = station_climate[tx_indexes, ])
        } else {
            tx_fit <- stats::lm(formula(paste0("tx", "~", "tn", "+", paste0(temps_covariates, collapse = "+"))),
                                data = station_climate[tx_indexes, ])
        }
        coefmax <- coefficients(tx_fit)

        if(control$use_stepwise) {
            tx_fit <- stats::step(tx_fit, trace = 0)
            coefmax[names(coefmax)] <- 0
            coefmax[names(coefficients(tx_fit))] <- coefficients(tx_fit)
        }

        station_climate[tx_indexes, "tx_residuals"] <- tx_fit$residuals

        # TODO: extract p-values for each covariate
        # coefs_summary <- summary(tx_fit)$coefficients
        # coefs_summary[, ncol(coefs_summary)]

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        tx_fit[["dates"]] <- station_climate[tx_indexes, "date"]
        # End of tx fit


        # Fit model for min temperature.
        tn_indexes <- na.omit(station_climate[, c("tn", "row_num", temps_covariates)])$row_num

        if (control$use_robust_methods) {
            tn_fit <- robustbase::lmrob(formula(paste0("tn", "~", paste0(temps_covariates, collapse = "+"))),
                                        data = station_climate[tn_indexes, ])
        } else {
            tn_fit <- stats::lm(formula(paste0("tn", "~", paste0(temps_covariates, collapse = "+"))),
                                data = station_climate[tn_indexes, ])
        }
        coefmin <- coefficients(tn_fit)

        if(control$use_stepwise) {
            tn_fit <- stats::step(tn_fit, trace = 0)
            coefmin[names(coefmin)] <- 0
            coefmin[names(coefficients(tn_fit))] <- coefficients(tn_fit)
        }

        station_climate[tn_indexes, "tn_residuals"] <- tn_fit$residuals

        # Se agregan las fechas para poder hacer las validaciones posteriormente
        tn_fit[["dates"]] <- station_climate[tn_indexes, "date"]
        # End of tn fit

        return_list <- list(coefficients = list(coefocc = coefocc, coefmin = coefmin, coefmax = coefmax),
                            gamma = list(coef = coefamt, alpha = alphamt),
                            residuals = station_climate[, c("date", "station", "probit_residuals", "tx_residuals", "tn_residuals")],
                            station = station_id)


        if(control$save_lm_fits) return_list[['lm_fits']] <- list(
            'tx' = tx_fit,
            'tn' = tn_fit,
            'prcp_occ' = prcp_occ_fit,
            'prcp_amount' = prcp_amt_fit
        )

        return(return_list)
    }

    first_model <- models[[1]]

    # Unwind results.  TODO: crear una func de .combine que haga esto de antemano y ahorrar este paso...
    models_coefficients <- list(
        coefocc = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefocc),
                         dimnames = list(gral.rownames, names(first_model$coefficients$coefocc))),

        coefmin = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmin),
                         dimnames = list(gral.rownames, names(first_model$coefficients$coefmin))),

        coefmax = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmax),
                         dimnames = list(gral.rownames, names(first_model$coefficients$coefmax))))

    models_residuals <- data.frame()

    GAMMA <- as.list(rep(NA, times = length(models)))
    names(GAMMA) <- gral.rownames

    for (m in models) {
        station <- as.character(m$station)
        models_coefficients$coefocc[station, ] <- m$coefficients$coefocc
        models_coefficients$coefmin[station, ] <- m$coefficients$coefmin
        models_coefficients$coefmax[station, ] <- m$coefficients$coefmax
        GAMMA[[station]] <- m$gamma
        models_residuals <- rbind(models_residuals, m$residuals)
    }

    model[['fit_metrics']] <- do.call(rbind, lapply(models, '[[', 'metrics'))
    model[["coefficients"]] <- models_coefficients
    model[["gamma"]] <- GAMMA

    if(control$save_lm_fits) {
        lm_fits <- lapply(models, '[[', 'lm_fits')
        names(lm_fits) <- unique_stations
        model[["lm_fits"]] <- lm_fits
    }

    rm(first_model, models, m)

    month_params <- foreach(m = 1:12, .multicombine = T, .export = c('partially_apply_LS'), .packages = c('dplyr')) %dopar% {
        month_residuals <- models_residuals %>% filter(lubridate::month(date) == m) %>% arrange(station, date)
        # month_residuals <- full_dates %>% left_join(month_residuals, by = 'date')
        month_climate <- climate %>% filter(lubridate::month(date) == m)

        # tx_residuals_matrix <- reshape2::acast(month_residuals, date ~ station, value.var = 'tx_residuals', fun.aggregate = function(x) ifelse(length(x) != 1, NA, x))
        tx_residuals_matrix <- reshape2::acast(month_residuals, date ~ station, value.var = 'tx_residuals')
        tn_residuals_matrix <- reshape2::acast(month_residuals, date ~ station, value.var = 'tn_residuals')
        probit_residuals_matrix <- reshape2::acast(month_residuals, date ~ station, value.var = 'probit_residuals')

        tx_matrix <- reshape2::acast(month_climate, date ~ station, value.var = 'tx')
        tn_matrix <- reshape2::acast(month_climate, date ~ station, value.var = 'tn')

        n_stations <- length(unique_stations)

        # tx_res_matrix <- matrix(month_residuals$tx_residuals, ncol = n_stations)
        # tn_res_matrix <- matrix(month_residuals$tn_residuals, ncol = n_stations)
        # probit_res_matrix <- matrix(month_residuals$probit_residuals, ncol = n_stations)

        # tx_matrix <- matrix(month_climate$tx, ncol = n_stations)
        # tn_matrix <- matrix(month_climate$tn, ncol = n_stations)

        # tx_matrix <- matrix(month_climate$tx, ncol=n_stations) tn_matrix <- matrix()

        # colnames(tx_residuals_matrix) <- colnames(tn_residuals_matrix) <- colnames(probit_residuals_matrix) <- colnames(tx_matrix) <- colnames(tn_matrix) <- unique_stations

        prcp_cor <- stats::cor(probit_residuals_matrix, use = control$cov_mode)

        prcp_vario <- stats::var(probit_residuals_matrix, use = control$cov_mode) * (1 - prcp_cor)
        tx_vario <- stats::cov(tx_residuals_matrix, use = control$cov_mode)
        tn_vario <- stats::cov(tn_residuals_matrix, use = control$cov_mode)

        # Assert that the column and row names of the variograms equal the ones of the distance matrix.
        stopifnot(all(colnames(prcp_vario) == gral.colnames), all(rownames(prcp_vario) == gral.rownames))
        stopifnot(all(colnames(tx_vario) == gral.colnames), all(rownames(tx_vario) == gral.rownames))
        stopifnot(all(colnames(tn_vario) == gral.colnames), all(rownames(tn_vario) == gral.rownames))

        if (fit_multiple_stations) {

            prcp_params <- stats::optimize(interval = c(0, max(dist.mat)), f = partially_apply_LS(prcp_vario, dist.mat, base_p = c(0, 1)))$minimum
            prcp_params <- c(0, 1, prcp_params)
            # prcp_params <- stats::optim(par = c(0.01, 1, max(dist.mat)), f = partially_apply_LS(prcp_vario, dist.mat))$par
            # prcp_params[params < 0] <- 0
            # prcp_params <- c(0, 1, prcp_params[3])

            # sill_initial_value <- mean(tx_vario[upper.tri(tx_vario, diag = T)])
            sill_initial_value <- mean(var(tx_matrix, na.rm = T))
            tx_params <- stats::optim(par = c(sill_initial_value, max(dist.mat)), fn = partially_apply_LS(tx_vario, dist.mat, base_p = c(0)))$par
            tx_params <- c(0, tx_params)
            # tx_params <- stats::optim(par = c(0.01, sill_initial_value, max(dist.mat)), fn = partially_apply_LS(tx_vario, dist.mat))$par
            # tx_params <- c(0, tx_params[2:3])

            # sill_initial_value <- mean(tn_vario[upper.tri(tn_vario, diag = T)])
            sill_initial_value <- mean(var(tn_matrix, na.rm = T))
            tn_params <- stats::optim(par = c(sill_initial_value, max(dist.mat)), fn = partially_apply_LS(tn_vario, dist.mat, base_p = c(0)))$par
            tn_params <- c(0, tn_params)

            # tn_params <- stats::optim(par = c(0.01, sill_initial_value, max(dist.mat)), fn = partially_apply_LS(tn_vario, dist.mat))$par
            # tn_params <- c(0, tn_params[2:3])

            variogram_parameters <- list(prcp = prcp_params, tx = tx_params, tn = tn_params)

        }

        residuals_sd <- data.frame(month_residuals %>% group_by(station) %>%
                                       summarise(tx = sd(tx_residuals, na.rm = T),
                                                 tn = sd(tn_residuals, na.rm = T),
                                                 prcp = 1))
        rownames(residuals_sd) <- residuals_sd[, 1]
        residuals_sd <- residuals_sd[, -1]

        cov_matrix <- list(tx = tx_vario, tn = tn_vario, prcp = prcp_cor)

        temp_ampl <- data.frame(month_climate %>% mutate(te = tx-tn) %>% group_by(station) %>%
                                    summarise(te_max = max(te, na.rm = T),
                                              te_min = min(te, na.rm = T)))
        rownames(temp_ampl) <- dplyr::pull(temp_ampl, station)
        temp_ampl <- dplyr::select(temp_ampl, -station)

        result <- list(residuals_sd = NULL,  cov_matrix = NULL, temp_ampl = NULL)
        if (fit_only_one_station)
            result <- list(residuals_sd = residuals_sd,
                           cov_matrix = cov_matrix,
                           temp_ampl = temp_ampl)
        if (fit_multiple_stations)
            result <- list(variogram_parameters = variogram_parameters,
                           residuals_sd = residuals_sd,
                           cov_matrix = cov_matrix,
                           temp_ampl = temp_ampl)

        return (result)
    }

    model[["month_params"]] <- month_params

    if(control$save_residuals) {
        model[["residuals"]] <- models_residuals
    }

    class(model) <- "glmwgen"

    model
}
