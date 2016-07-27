partial <- function(f, ...) {
    l <- list(...)
    function(...) {
        do.call(f, c(list(...), l))
    }
}

get_seasonal_covariat <- function(season_number, season_values, break_quantiles = c(0, 0.33, 0.66, 1), levels_probabilities = c(1/3, 1/3, 1/3)) {
    splitted <- cut(season_values, breaks = quantile(season_values, probs = break_quantiles), include.lowest = T)
    values_level <- sample(levels(splitted), 1, prob = levels_probabilities)
    return(sample(season_values[splitted == values_level], 1))
}

get_rainfall_seasonal_covariat <- get_seasonal_covariat
get_temperatures_seasonal_covariat <- function(season_number, season_tx_values, season_tn_values, break_quantiles = c(0, 0.33, 0.66, 1),
                                               levels_probabilities = c(0.25, 0.35, 0.4)) {
    splitted_tx <- cut(season_tx_values, breaks = quantile(season_tx_values, probs = break_quantiles), include.lowest = T)
    splitted_tn <- cut(season_tn_values, breaks = quantile(season_tn_values, probs = break_quantiles), include.lowest = T)

    values_level_index <- sample(seq_along(levels(splitted_tx)), 1, prob = levels_probabilities)
    tx_sampled_level <- levels(splitted_tx)[values_level_index]
    tn_sampled_level <- levels(splitted_tn)[values_level_index]
    return(list(tx = sample(season_tx_values[splitted_tx == tx_sampled_level], 1),
                tn = sample(season_tn_values[splitted_tn == tn_sampled_level], 1)))
}

#' @title Simulations control configuration
#' @description Provides fine control of different parameters that will be used to create new weather series.
#' @export
glmwgen_simulation_control <- function(seasonal_temps_covariats_getter = get_temperatures_seasonal_covariat,
                                     seasonal_prcp_covariats_getter = get_rainfall_seasonal_covariat,
                                     random_fields_method = 'gaussian',
                                     Rt = NULL,
                                     grf_method = NULL,
                                     multicombine = T) {
    rf_function <- NULL
    if(!is.function(random_fields_method)) {
        if(startsWith(random_fields_method, 'gauss')) rf_function <- gaussian_random_field
        if(startsWith(random_fields_method, 'chol')) rf_function <- cholesky_random_field
    } else {
        rf_function <- random_fields_method
    }
    return(list(seasonal_temps_covariats_getter = seasonal_temps_covariats_getter,
                seasonal_prcp_covariats_getter = seasonal_prcp_covariats_getter,
                random_fields_method = rf_function,
                Rt = Rt,
                grf_method = grf_method,
                multicombine = multicombine))
}

#' @title Simulates new weather trajectories
#' @description Simulates new weather trajectories.
#' @param object A glmwgen model.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#'
#'          Either NULL or an integer that will be used in a call to set.seed before simulating the response vectors.
#'          If set, the value is saved as the "seed" attribute of the returned value. The default, NULL will not change the random generator state.
#' @param start_date a start date in text format (will be converted using as.Date) or a date object.
#' @param end_date an end date in text format (will be converted using as.Date) or a date object.
#' @param simulation_locations a SpatialPointsDataFrame with the points at which weather should be simulated.
#'          If not set, the locations used to fit the model will be used.
#' @param control a glmwgen simulation control list.
#' @import dplyr
#' @import foreach
#' @export
simulate.glmwgen <- function(object, nsim = 1, seed = NULL, start_date = NA, end_date = NA, simulation_locations = NULL, control = glmwgen_simulation_control(),
                             verbose = T) {
    model <- object

    if(class(object) != 'glmwgen') {
        stop(paste('Received a model of class', class(object), 'and a model of class "glmwgen" was expected.'))
    }

    if(!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(simulation_locations)) {
        simulation_locations <- coordinates(model$stations)
    }

    simulation_dates <- data.frame(date = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "days"))

    if (is.null(control$Rt)) {
        # TODO: shouldn't this be = 1 once we start simulating?
        Rt <- seq(from = -1, to = 1, length.out = length(simulation_dates))
    } else {
        Rt <- control$Rt
        if (length(Rt) == 1)
            Rt <- rep(Rt, times = length(simulation_dates))
        if (length(Rt) != length(simulation_dates))
            stop("The specified Rt parameter length differs from the simulation dates series length.")
    }

    simulation_dates <- simulation_dates %>%
        mutate(year_fraction = 2 * pi * lubridate::yday(date) / ifelse(lubridate::leap_year(date), 366, 365),
               ct = cos(year_fraction),
               st = sin(year_fraction),
               Rt = Rt,
               season = ceiling(lubridate::month(date)/3),
               ST1 = 0, ST2 = 0, ST3 = 0, ST4 = 0,
               SMX1 = 0, SMX2 = 0, SMX3 = 0, SMX4 = 0,
               SMN1 = 0, SMN2 = 0, SMN3 = 0, SMN4 = 0)

    rm(Rt)

    # estimate model coefficients on grid using ordinary kriging (OK) occurrence
    coefocc.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefocc),
                          dimnames = list(NULL, colnames(model$coefficients$coefocc)))
    for (covariat in 1:ncol(coefocc.sim)) {
        coefocc.sim[, covariat] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                         model$coefficients$coefocc[, covariat]),
                                                            simulation_locations))
    }

    # minimum temperature
    coefmin.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefmin),
                          dimnames = list(NULL, colnames(model$coefficients$coefmin)))
    for (covariat in 1:ncol(coefmin.sim)) {
        coefmin.sim[, covariat] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                         model$coefficients$coefmin[, covariat]),
                                                            simulation_locations))
    }

    # maximum temperature
    coefmax.sim <- matrix(NA, nrow = nrow(simulation_locations), ncol = ncol(model$coefficients$coefmax),
                          dimnames = list(NULL, colnames(model$coefficients$coefmax)))
    for (covariat in 1:ncol(coefmax.sim)) {
        coefmax.sim[, covariat] <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations),
                                                                         model$coefficients$coefmax[, covariat]),
                                                            simulation_locations))
    }

    # Gamma shape and scale
    SH <- c()

    for (station in as.character(model$stations$id)) SH[station] <- model$gamma[[station]]$alpha
    SH.sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SH), simulation_locations))
    # SH.sim <- fields::predict.Krig(fields::Krig(coordinates(coordinates(model$stations)), SH), grid_test)
    if(verbose) cat('Kriged coefficients..\n')


    if(identical(control$random_fields_method, cholesky_random_field)) {
        cat('Calculate cSigma matrixes.\n')
        # cSigma_list <- list()
        # # grid     <- list(x = unique(coordinates(simulation_locations)[,1]), y = unique(coordinates(simulation_locations)[,2]))
        # # xg       <- fields::make.surface.grid(grid)
        #
        # model$cSigma <- foreach(month_number=1:12, .multicombine = T) %do% {
        #     month_params <- model$month_params[[month_number]]
        #     cSigma_list <- list()
        #     for(var_name in names(month_params)) {
        #
        #         # bigSigma <- fields::stationary.cov(xg, theta=month_params[[var_name]][3])
        #         bigSigma <- fields::stationary.cov(simulation_locations, theta=month_params[[var_name]][3])
        #         cSigma_list[[var_name]] <- t(chol(bigSigma))
        #     }
        #     print(paste0('Chol transform. Month = ', month, '. Variable = ', var_name))
        #
        #     cSigma_list
        # }
        #
        # # names(model$cSigma) <- unique(lubridate::month(simulation_dates$date))
        #
        # # rm(grid, xg, bigSigma)
        # rm(bigSigma)
    }


    # gen_climate <- foreach(i=1:nsim, .combine=partial(abind, along=0.5), .multicombine=T) %dopar% {
    gen_climate <- foreach(i = 1:nsim, .combine = list, .multicombine = control$multicombine) %do% {
        simulation_start <- as.Date(start_date) - 1

        # initialize with climatology of the previous day
        # w2 <- suppressWarnings(geoR::grf(nrow(simulation_locations),
        #                                  grid = simulation_locations,
        #                                  cov.model = "exponential",
        #                                  cov.pars = model$month_params[[lubridate::month(simulation_start)]]$prcp[c(2,3)],
        #                                  nugget = model$month_params[[lubridate::month(simulation_start)]]$prcp[1],
        #                                  mean = rep(0, nrow(simulation_locations)), messages = FALSE))
        w2 <- control$random_fields_method(model, simulation_locations, control, lubridate::month(simulation_start), 'prcp')
        SIMocc.old <- (w2 > 0) + 0

        start_climatology <- model$start_climatology %>% filter(month == lubridate::month(simulation_start), day == lubridate::day(simulation_start))

        SIMmin.old <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tn), simulation_locations)))
        SIMmax.old <- suppressWarnings(as.vector(fields::predict.Krig(fields::Krig(coordinates(model$stations), start_climatology$tx), simulation_locations)))

        SC <- array(data = NA, dim = c(nrow(simulation_dates), nrow(coordinates(model$stations))))
        colnames(SC) <- model$stations$id

        for (season_number in unique(simulation_dates$season)) {
            season_indexes <- simulation_dates$season == season_number
            # season_values <- d_seatot[d_seatot$season == season_number, ]
            season_tx <- unique(model$seasonal$tx[[season_number]])
            season_tn <- unique(model$seasonal$tn[[season_number]])
            season_prcp <- unique(model$seasonal$prcp[[season_number]])
            temp_covariats <- control$seasonal_temps_covariats_getter(season_indexes, season_tx[season_tx > 0], season_tn[season_tn > 0])

            simulation_dates[season_indexes, paste0("ST", season_number)] <- control$seasonal_prcp_covariats_getter(season_number, season_prcp[season_prcp > 0])
            simulation_dates[season_indexes, paste0("SMX", season_number)] <- temp_covariats$tx
            simulation_dates[season_indexes, paste0("SMN", season_number)] <- temp_covariats$tn
        }

        daily_covariats <- as.matrix(t(simulation_dates[, c(3:5, 7:18)]))

        daily_covariats <- rbind(1, daily_covariats)
        rownames(daily_covariats)[1] <- "(Intercept)"

        for (station in as.character(model$stations$id)) {
            gamma_coef <- model$gamma[[station]]$coef
            SC[, station] <- exp(apply(daily_covariats[names(gamma_coef), ] * gamma_coef, FUN = sum, MAR = 2, na.rm = T))/SH[station]
        }

        simulated_occurrence <- simulated_tx <- simulated_tn <- simulated_prcp <- array(dim = c(nrow(simulation_dates), nrow(simulation_locations)))
        # dimnames = list(c('date', 'location')))

        if(verbose) cat('Simulated start climatology.\n')

        temps_retries <- 0
        for (d in 1:nrow(simulation_dates)) {
            # month_params <- model$month_params[[lubridate::month(simulation_start)]]
            month_number <- lubridate::month(simulation_dates[d, "date"])
            # p <- model$month_params[[lubridate::month(simulation_dates[d, "date"])]]

            simulation_matrix <- cbind(SIMocc.old, SIMmin.old, SIMmax.old, matrix(daily_covariats[, d], ncol = nrow(daily_covariats), nrow = length(SIMmax.old), byrow = T))
            colnames(simulation_matrix) <- c("prcp_occ_prev", "tn_prev", "tx_prev", rownames(daily_covariats))


            mu.occ <- apply(simulation_matrix[, colnames(coefocc.sim)] * coefocc.sim, 1, sum, na.rm = T)
            # w.occ <- suppressWarnings(geoR::grf(nrow(simulation_locations), grid = simulation_locations, cov.model = "exponential", cov.pars = c(p$prcp[2], p$prcp[3]), nugget = p$prcp[1],
            #     mean = rep(0, nrow(simulation_locations)), messages = FALSE)$data)
            w.occ <- control$random_fields_method(model, simulation_locations, control, month_number, 'prcp')

            simulated_occurrence[d, ] <- ((mu.occ + w.occ) > 0) + 0

            simulation_matrix <- cbind(prcp_occ = simulated_occurrence[d, ], simulation_matrix)

            mu.min <- apply(simulation_matrix[, colnames(coefmin.sim)] * coefmin.sim, 1, sum, na.rm = T)
            mu.max <- apply(simulation_matrix[, colnames(coefmax.sim)] * coefmax.sim, 1, sum, na.rm = T)

            # w.min <- suppressWarnings(geoR::grf(nrow(simulation_locations),
            #                                     grid = simulation_locations,
            #                                     cov.model = "exponential",
            #                                     cov.pars = c(p$tn[2], p$tn[3]),
            #                                     nugget = p$tn[1],
            #                                     mean = rep(0, nrow(simulation_locations)),
            #                                     messages = FALSE)$data)
            w.min <- control$random_fields_method(model, simulation_locations, control, month_number, 'tn')
            w.max <- control$random_fields_method(model, simulation_locations, control, month_number, 'tx')

            # w.max <- suppressWarnings(geoR::grf(nrow(simulation_locations),
            #                                     grid = simulation_locations,
            #                                     cov.model = "exponential",
            #                                     cov.pars = c(p$tx[2], p$tx[3]),
            #                                     nugget = p$tx[1],
            #                                     mean = rep(0, nrow(simulation_locations)),
            #                                     messages = FALSE)$data)

            simulated_tn[d, ] <- signif(mu.min + w.min, digits = 4)
            simulated_tx[d, ] <- signif(mu.max + w.max, digits = 4)

            # 1/10000
            while (min(simulated_tx[d, ] - simulated_tn[d, ]) < 0.5) {
                # if(verbose) print(paste('Retrying min and max temperatures simulation.', 'Min difference: ', min(simulated_tx[d, ] - simulated_tn[d, ])))
                temps_retries <- temps_retries + 1

                # w.min <- suppressWarnings(geoR::grf(nrow(simulation_locations), grid = simulation_locations, cov.model = "exponential", cov.pars = c(p$tn[2], p$tn[3]), nugget = p$tn[1],
                #   mean = rep(0, nrow(simulation_locations)), messages = FALSE)$data)
                #
                #
                # w.max <- suppressWarnings(geoR::grf(nrow(simulation_locations), grid = simulation_locations, cov.model = "exponential", cov.pars = c(p$tx[2], p$tx[3]), nugget = p$tx[1],
                #   mean = rep(0, nrow(simulation_locations)), messages = FALSE)$data)

                w.min <- control$random_fields_method(model, simulation_locations, control, month_number, 'tn')
                w.max <- control$random_fields_method(model, simulation_locations, control, month_number, 'tx')

                simulated_tn[d, ] <- signif(mu.min + w.min, digits = 4)
                simulated_tx[d, ] <- signif(mu.max + w.max, digits = 4)
            }

            SIMocc.old <- simulated_occurrence[d, ]
            SIMmax.old <- simulated_tx[d, ]
            SIMmin.old <- simulated_tn[d, ]

            ## amounts
            SC.sim <- suppressWarnings(fields::predict.Krig(fields::Krig(coordinates(model$stations), SC[d, ]), simulation_locations))

            # w3 <- suppressWarnings(geoR::grf(nrow(simulation_locations), grid = simulation_locations, cov.model = "exponential", cov.pars = c(p$prcp[2], p$prcp[3]), nugget = p$prcp[1], mean = rep(0,
            #     nrow(simulation_locations)), messages = FALSE)$data)
            w3 <- control$random_fields_method(model, simulation_locations, control, month_number, 'prcp')

            simulated_prcp[d, ] <- signif(qgamma(pnorm(w3), shape = SH.sim, scale = SC.sim), digits = 4)
            if(verbose) cat(paste0("\r", nsim, ": ", d, "/", ncol(daily_covariats), ". Retries: ", temps_retries))
        }

        simulated_prcp[simulated_prcp < model$control$prcp_occurrence_threshold] <- 0
        simulated_occurrence[simulated_prcp == 0] <- 0
        simulated_prcp[simulated_occurrence == 0] <- 0

        return(list(prcp = simulated_prcp, tx = simulated_tx, tn = simulated_tn))
    }

    # save(gen_climate, file='gen_climate.RData')


    if (nsim == 1) {
        gen_climate <- list(gen_climate)
    }

    gen_climate
}
