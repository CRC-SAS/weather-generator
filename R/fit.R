#' @title Weather model fit configuration
#' @description Provides fine control of different parameters that will be used to fit a weather model.
#' @export
glmwgenControl <- function(prcp_occurrence_threshold = 0.1, use_seasonal_covariats = F, seasonal_covariats = c("tx", "tn", "prcp")) {
    if (!all(seasonal_covariats %in% c("tx", "tn", "prcp"))) {
        stop("The seasonal_covariats parameter should list the variables names that will be fitted with seasonal averages as covariats.")
    }
    return(list(prcp_occurrence_threshold = prcp_occurrence_threshold,
                use_seasonal_covariats = use_seasonal_covariats,
                seasonal_covariats = seasonal_covariats))
}


#' @title Weather model calibration
#' @description Fits a weather model from historical data.
#' @import foreach
#' @import dplyr
#' @export
calibrate.glmwgen <- function(climate, stations, control = glmwgenControl(), verbose = T) {
    model <- list()
    climate <- climate %>% arrange(station, date)
    stations <- stations[order(stations$id), ]

    unique_stations <- unique(climate$station)

    if (!all(unique_stations == stations$id)) {
        stop("Mismatch between stations ids between climate and stations datasets.")
    }

    # distance matrix of station locations converts longitude and latitude to KILOMETERS
    # dist.mat <- fields::rdist.earth(coordinates(stations), miles = F)
    dist.mat <- sp::spDists(stations)
    # dist.mat <- sp::spDists(spTransform(stations, CRS("+proj=longlat +datum=WGS84")))
    colnames(dist.mat) <- stations$id
    rownames(dist.mat) <- stations$id
    # diagonal element is not explicitly equal to zero, so define as such
    diag(dist.mat) <- 0

    seasonal_covariats <- estimate_seasonal_covariats(climate)

    ## TODO: check control variables.
    model[["control"]] <- control

    model[["distance_matrix"]] <- dist.mat
    model[["seasonal"]] <- seasonal_covariats
    model[['stations_proj4string']] <- stations@proj4string
    model[["stations"]] <- stations
    # model[["stations"]] <- sp::spTransform(stations, CRS("+proj=longlat +datum=WGS84"))

    model[["start_climatology"]] <- climate %>%
        group_by(station, month = lubridate::month(date), day = lubridate::day(date)) %>%
        summarise(tx = mean(tx, na.rm = T),
                  tn = mean(tn, na.rm = T),
                  prcp = mean(prcp, na.rm = T))

    prcp_covariats <- c("ct", "st")
    temps_covariats <- c("tn_prev", "tx_prev", "ct", "st", "prcp_occ", "Rt")

    if (control$use_seasonal_covariats) {
        prcp_covariats <- c(prcp_covariats, "ST1", "ST2", "ST3", "ST4")
        temps_covariats <- c(temps_covariats, "SMN1", "SMN2", "SMN3", "SMN4", "SMX1", "SMX2", "SMX3", "SMX4")
    }

    models <- foreach::foreach(station = unique_stations, .multicombine = T) %dopar% {
        # system.time(replicate(500, {
        station_climate <- climate[climate$station == station, ] %>%
            mutate(prcp_occ = (prcp > control$prcp_occurrence_threshold) + 0,
                   prcp_intensity = ifelse(prcp < control$prcp_occurrence_threshold, NA, prcp),
                   prcp_occ_prev = lag(prcp_occ),
                   prcp_prev = lag(prcp),
                   tx_prev = lag(tx),
                   tn_prev = lag(tn),
                   year_fraction = 2 * pi * lubridate::yday(date)/ifelse(lubridate::leap_year(date), 366, 365),
                   ct = cos(year_fraction),
                   st = sin(year_fraction),
                   Rt = seq(from = -1, to = 1, length.out = length(year_fraction)),
                   row_num = row_number())



        if (control$use_seasonal_covariats) {
            station_climate <- data.frame(station_climate,
                                          ST1 = seasonal_covariats$prcp[[1]],
                                          ST2 = seasonal_covariats$prcp[[2]],
                                          ST3 = seasonal_covariats$prcp[[3]],
                                          ST4 = seasonal_covariats$prcp[[4]],
                                          SMX1 = seasonal_covariats$tx[[1]],
                                          SMX2 = seasonal_covariats$tx[[2]],
                                          SMX3 = seasonal_covariats$tx[[3]],
                                          SMX4 = seasonal_covariats$tx[[4]],
                                          SMN1 = seasonal_covariats$tn[[1]],
                                          SMN2 = seasonal_covariats$tn[[2]],
                                          SMN3 = seasonal_covariats$tn[[3]],
                                          SMN4 = seasonal_covariats$tn[[4]])
        }

        # Fit model for precipitation occurrence.
        probit_indexes <- na.omit(station_climate[, c("prcp_occ", "row_num", "prcp_occ_prev", prcp_covariats)])$row_num

        occ_fit <- stats::glm(formula(paste0("prcp_occ", "~", "prcp_occ_prev +", paste0(prcp_covariats, collapse = "+"))),
                             data = station_climate[probit_indexes, ],
                             family = stats::binomial(probit))

        coefocc <- occ_fit$coefficients
        station_climate[probit_indexes, "probit_residuals"] <- occ_fit$residuals

        # Fit model for precipitation amounts.
        gamma_indexes <- na.omit(station_climate[, c("prcp_intensity", "row_num", prcp_covariats)])$row_num
        # save model in list element 'i'
        prcp_fit <- stats::glm(formula(paste0("prcp_intensity", "~", paste0(prcp_covariats, collapse = "+"))),
                              data = station_climate[gamma_indexes, ],
                              family = stats::Gamma(link = log))

        # coefamt <- prcp_fit$coefficients

        tx_indexes <- na.omit(station_climate[, c("tx", "row_num", temps_covariats)])$row_num
        tn_indexes <- na.omit(station_climate[, c("tn", "row_num", temps_covariats)])$row_num
        # Fit

        tx_fit <- stats::lm(formula(paste0("tx", "~", paste0(temps_covariats, collapse = "+"))), data = station_climate[tx_indexes, ])
        coefmax <- tx_fit$coefficients
        station_climate[tx_indexes, "tx_residuals"] <- tx_fit$residuals

        tn_fit <- stats::lm(formula(paste0("tn", "~", paste0(temps_covariats, collapse = "+"))), data = station_climate[tn_indexes, ])
        coefmin <- tn_fit$coefficients
        station_climate[tn_indexes, "tn_residuals"] <- tn_fit$residuals

        return(list(coefficients = list(coefocc = coefocc, coefmin = coefmin, coefmax = coefmax),
                    gamma = list(coef = prcp_fit$coef, alpha = MASS::gamma.shape(prcp_fit)$alpha),
                    residuals = station_climate[, c("date", "station", "probit_residuals", "tx_residuals", "tn_residuals")],
                    station = station))

        # }))
    }

    first_model <- models[[1]]

    # Unwind results.  TODO: crear una func de .combine que haga esto de antemano y ahorrar este paso...
    models_coefficients <- list(
        coefocc = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefocc),
                         dimnames = list(rownames(dist.mat), names(first_model$coefficients$coefocc))),

        coefmin = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmin),
                         dimnames = list(rownames(dist.mat), names(first_model$coefficients$coefmin))),

        coefmax = matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmax),
                         dimnames = list(rownames(dist.mat), names(first_model$coefficients$coefmax))))

    models_residuals <- data.frame()

    GAMMA <- as.list(rep(NA, times = length(models)))
    names(GAMMA) <- rownames(dist.mat)

    for (m in models) {
        station <- as.character(m$station)
        models_coefficients$coefocc[station, ] <- m$coefficients$coefocc
        models_coefficients$coefmin[station, ] <- m$coefficients$coefmin
        models_coefficients$coefmax[station, ] <- m$coefficients$coefmax
        GAMMA[[station]] <- m$gamma
        models_residuals <- rbind(models_residuals, m$residuals)
    }

    model[["coefficients"]] <- models_coefficients
    model[["gamma"]] <- GAMMA

    rm(first_model, models, m)

    month_params <- foreach(m = 1:12, .multicombine = T) %dopar% {
        month_residuals <- models_residuals %>% filter(lubridate::month(date) == m) %>% arrange(station, date)
        month_climate <- climate %>% filter(lubridate::month(date) == m)

        n_stations <- length(unique_stations)

        tx_res_matrix <- matrix(month_residuals$tx_residuals, ncol = n_stations)
        tn_res_matrix <- matrix(month_residuals$tn_residuals, ncol = n_stations)
        probit_res_matrix <- matrix(month_residuals$probit_residuals, ncol = n_stations)

        # tx_matrix <- matrix(month_climate$tx, ncol=n_stations) tn_matrix <- matrix()

        colnames(tx_res_matrix) <- colnames(tn_res_matrix) <- colnames(probit_res_matrix) <- unique_stations


        # all(na.omit(tx_matrix[, '87448']) == na.omit(month_climate[month_climate$station == 87448, 'tx']))
        # all(na.omit(tx_res_matrix[,1]) == na.omit(month_residuals[month_residuals$station == 87448, 'tx_residuals']))

        prcp_cor <- stats::cor(probit_res_matrix, use = "pairwise.complete")

        PRCPvario <- stats::var(probit_res_matrix, use = "pairwise.complete") * (1 - prcp_cor)
        TMAXvario <- stats::cov(tx_res_matrix, use = "pairwise.complete")
        TMINvario <- stats::cov(tn_res_matrix, use = "pairwise.complete")

        # Assert that the column and row names of the variograms equal the ones of the distance matrix.
        stopifnot(all(colnames(TMAXvario) == colnames(dist.mat)), all(rownames(TMAXvario) == rownames(dist.mat)))

        params <- stats::optimize(interval = c(0, max(dist.mat)), f = partially_apply_LS(PRCPvario, dist.mat, base_p = c(0, 1)))$minimum
        params <- c(0, 1, params)

        sill_initial_value <- mean(TMAXvario[upper.tri(TMAXvario, diag = T)])
        params.max <- stats::optim(par = c(sill_initial_value, max(dist.mat)), fn = partially_apply_LS(TMAXvario, dist.mat, base_p = c(0)))$par
        params.max <- c(0, params.max)

        sill_initial_value <- mean(TMINvario[upper.tri(TMINvario, diag = T)])
        params.min <- stats::optim(par = c(sill_initial_value, max(dist.mat)), fn = partially_apply_LS(TMINvario, dist.mat, base_p = c(0)))$par
        params.min <- c(0, params.min)

        return(list(prcp = params, tx = params.max, tn = params.min))
    }

    model[["month_params"]] <- month_params

    rm(models_residuals)

    class(model) <- "glmwgen"

    model
}
