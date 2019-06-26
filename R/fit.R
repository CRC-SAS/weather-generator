
#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
glmwgen_fit_control <- function(prcp_occurrence_threshold = 0.1,
                           use_seasonal_covariates = F, # el modelo más simple es sin covariables
                           use_linear_term = T,
                           seasonal_covariates = c("tx", "tn", "prcp"),
                           save_residuals = F,
                           cov_mode = 'complete',
                           save_lm_fits = F,
                           use_stepwise = F,
                           use_robust_methods = F) {

    if (!all(seasonal_covariates %in% c("tx", "tn", "prcp"))) {
        stop("The seasonal_covariates parameter should list the variables names that will be fitted with seasonal averages as covariates.")
    }

    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                use_seasonal_covariates = use_seasonal_covariates,
                use_linear_term = use_linear_term,
                seasonal_covariates = seasonal_covariates,
                save_residuals = save_residuals,
                cov_mode = cov_mode,
                save_lm_fits = save_lm_fits,
                use_stepwise = use_stepwise,
                use_robust_methods = use_robust_methods))
}


#' @title Weather model calibration
#' @description Fits a weather model from historical data.
#' @import foreach
#' @import dplyr
#' @export
calibrate.glmwgen <- function(climate, stations, control = glmwgen_fit_control(), verbose = T) {

    check.input.data(climate, stations)

    model <- list()

    climate <- climate %>%
        dplyr::rename(station = station_id)
    climate <- climate %>% arrange(station, date)

    stations <- stations %>%
        dplyr::rename(id = station_id) %>%
        sf::as_Spatial()
    stations <- stations[order(stations$id), ]
    rownames(stations@data) <- stations$id

    fit_multiple_stations <- fit_only_one_station <- F
    if (nrow(stations) == 1) fit_only_one_station <- T
    if (nrow(stations) > 1) fit_multiple_stations <- T

    gral.rownames <- as.character(stations$id)
    gral.colnames <- as.character(stations$id)

    unique_stations <- sort(unique(climate$station))

    if (!all(unique_stations == stations$id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    invalid_records <- c(which(climate$tx < climate$tn), which(climate$prcp < 0))

    if(length(invalid_records) > 0) {
        warning(sprintf('Removing %d records with maximum temperature < minimum temperature or negative rainfalls.', length(invalid_records)))
        climate[invalid_records, c('tx', 'tn', 'prcp')] <- NA
    }

    if (fit_multiple_stations) {
        dist.mat <- sp::spDists(stations) # no se debe usar al argumento longlat porque fit se hace diferentes proyecciones!
        colnames(dist.mat) <- gral.colnames
        rownames(dist.mat) <- gral.rownames
        # diagonal element is not explicitly equal to zero, so define as such
        diag(dist.mat) <- 0
    }

    summarised_climate <- summarise_seasonal_climate(climate)

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

    prcp_covariates <- c("ct", "st")
    temps_covariates <- c("tn_prev", "tx_prev", "ct", "st", "prcp_occ")

    if(control$use_linear_term) {
        temps_covariates <- c(temps_covariates, "Rt")
    }

    if (control$use_seasonal_covariates) {
        prcp_covariates <- c(prcp_covariates, "ST1", "ST2", "ST3", "ST4")
        temps_covariates <- c(temps_covariates, "SMN1", "SMN2", "SMN3", "SMN4", "SMX1", "SMX2", "SMX3", "SMX4")
    }

    # Register a sequential backend if the user didn't register a parallel
    # in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        foreach::registerDoSEQ()
    }

    full_dates <- data.frame(date = seq.Date(min(climate$date), max(climate$date), by = 'days'))

    models <- foreach::foreach(station = unique_stations, .multicombine = T, .packages = c('dplyr')) %dopar% {
        station_climate <- full_dates %>% left_join(climate[climate$station == station, ], by = 'date') %>%
            mutate(prcp_occ = (prcp >= control$prcp_occurrence_threshold) + 0,              # prcp occurrence
                   prcp_amt = ifelse(prcp >= control$prcp_occurrence_threshold, prcp, NA),  # prcp amount/intensity
                   prcp_occ_prev = lag(prcp_occ),
                   tx_prev = lag(tx),
                   tn_prev = lag(tn),
                   year_fraction = 2 * pi * lubridate::yday(date)/ifelse(lubridate::leap_year(date), 366, 365),
                   ct = cos(year_fraction),
                   st = sin(year_fraction),
                   Rt = seq(from = -1, to = 1, length.out = length(year_fraction)),
                   row_num = row_number())

        station_climate$station <- station



        if (control$use_seasonal_covariates) {
            seasonal_covariates <- create_seasonal_covariates(summarised_climate)

            station_climate <- data.frame(station_climate,
                                          ST1 = seasonal_covariates$prcp[[1]],
                                          ST2 = seasonal_covariates$prcp[[2]],
                                          ST3 = seasonal_covariates$prcp[[3]],
                                          ST4 = seasonal_covariates$prcp[[4]],
                                          SMX1 = seasonal_covariates$tx[[1]],
                                          SMX2 = seasonal_covariates$tx[[2]],
                                          SMX3 = seasonal_covariates$tx[[3]],
                                          SMX4 = seasonal_covariates$tx[[4]],
                                          SMN1 = seasonal_covariates$tn[[1]],
                                          SMN2 = seasonal_covariates$tn[[2]],
                                          SMN3 = seasonal_covariates$tn[[3]],
                                          SMN4 = seasonal_covariates$tn[[4]])
        }

        # Fit model for precipitation occurrence.
        probit_indexes <- na.omit(station_climate[, c("prcp_occ", "row_num", "prcp_occ_prev", prcp_covariates)])$row_num

        if (control$use_robust_methods) {
            prcp_occ_fit <- robust::glmRob(formula(paste0("prcp_occ", "~", "prcp_occ_prev +", paste0(prcp_covariates, collapse = "+"))),
                                      data = station_climate[probit_indexes, ],
                                      family = stats::binomial(probit))
        } else {
            prcp_occ_fit <- stats::glm(formula(paste0("prcp_occ", "~", "prcp_occ_prev +", paste0(prcp_covariates, collapse = "+"))),
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

        prcp_occ_fit_auc <- NA
        if("ROCR" %in% rownames(installed.packages())) {
            tmp_auc <- ROCR::prediction(fitted(prcp_occ_fit), (na.omit(station_climate[, c("prcp_occ", "prcp_occ_prev", prcp_covariates)]))$prcp_occ)
            prcp_occ_fit_auc <- ROCR::performance(tmp_auc, 'auc')@y.values[[1]]
            rm(tmp_auc)
        } else {
            warning ("Install ROCR package to get AUC metrics for the precipitation occurrences fit!")
        }
        prcp_occ_fit_rsq <- NA
        if("ResourceSelection" %in% rownames(installed.packages())) {
            tmp_rsq <- ResourceSelection::hoslem.test(prcp_occ_fit$y, prcp_occ_fit$fitted.values, g = 100)
            prcp_occ_fit_rsq <- tmp_rsq$statistic
            rm(tmp_rsq)
        } else {
            warning ("Install ResourceSelection package to get R² metrics for the precipitation occurrences fit!")
        }


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

        tx_indexes <- na.omit(station_climate[, c("tx", "row_num", temps_covariates)])$row_num
        tn_indexes <- na.omit(station_climate[, c("tn", "row_num", temps_covariates)])$row_num
        # Fit

        if (control$use_robust_methods) {
            tx_fit <- robustbase::lmrob(formula(paste0("tx", "~", paste0(temps_covariates, collapse = "+"))),
                                        data = station_climate[tx_indexes, ])
        } else {
            tx_fit <- stats::lm(formula(paste0("tx", "~", paste0(temps_covariates, collapse = "+"))),
                                data = station_climate[tx_indexes, ])
        }
        coefmax <- coefficients(tx_fit)

        if(control$use_stepwise) {
            tx_fit <- stats::step(tx_fit, trace = 0)
            coefmax[names(coefmax)] <- 0
            coefmax[names(coefficients(tx_fit))] <- coefficients(tx_fit)
        }

        station_climate[tx_indexes, "tx_residuals"] <- tx_fit$residuals

        tx_fit_rsq <- summary(tx_fit)$r.squared

        # TODO: extract p-values for each covariate
        # coefs_summary <- summary(tx_fit)$coefficients
        # coefs_summary[, ncol(coefs_summary)]

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

        tn_fit_rsq <- summary(tn_fit)$r.squared

        return_list <- list(coefficients = list(coefocc = coefocc, coefmin = coefmin, coefmax = coefmax),
                            gamma = list(coef = coefamt, alpha = alphamt),
                            residuals = station_climate[, c("date", "station", "probit_residuals", "tx_residuals", "tn_residuals")],
                            metrics = data.frame(station = station,
                                                 metric = c('R²', 'R²', 'AUC', 'R²'),
                                                 model = c('tx', 'tn', 'prcp_occ', 'prcp_occ'),
                                                 value = c(tx_fit_rsq, tn_fit_rsq, prcp_occ_fit_auc, prcp_occ_fit_rsq)),
                            # significance = list(coefocc = )
                            station = station)


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

        result <- list(residuals_sd = NULL,  cov_matrix = NULL)
        if (fit_only_one_station)
            result <- list(residuals_sd = residuals_sd,  cov_matrix = cov_matrix)
        if (fit_multiple_stations)
            result <- list(variogram_parameters = variogram_parameters, residuals_sd = residuals_sd,  cov_matrix = cov_matrix)

        return (result)
    }

    model[["month_params"]] <- month_params

    if(control$save_residuals) {
        model[["residuals"]] <- models_residuals
    }

    class(model) <- "glmwgen"

    model
}
