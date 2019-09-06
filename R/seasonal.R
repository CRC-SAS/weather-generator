
summarise_seasonal_climate <- function(datos_climaticos, umbral.faltantes = 0.2) {
    ## compute seasonal averages for tx and tn, and sum for prcp, with quality control and missing values imputation
    purrr::map_dfr(
        .x = unique(dplyr::pull(datos_climaticos, station)),
        .f = function(omm.id) {
            datos.omm.id     <- dplyr::filter(datos_climaticos, station == omm.id) %>%
                tidyr::complete(date = base::seq.Date(from = lubridate::floor_date(min(date), "year") ,
                                                       to   = lubridate::ceiling_date(max(date), "year") - lubridate::days(1),
                                                       by   = "days")) %>%
                dplyr::mutate(season = lubridate::quarter(date, fiscal_start = 12), year = lubridate::year(date))
            estadisticas     <- datos.omm.id %>%
                dplyr::group_by(station, year, season) %>%
                dplyr::summarise(cantidad_datos = dplyr::n(),
                                 cantidad_tmax = sum(ifelse(is.na(tx), 0, 1)),
                                 cantidad_tmin = sum(ifelse(is.na(tn), 0, 1)),
                                 cantidad_prcp = sum(ifelse(is.na(prcp), 0, 1)),
                                 mean_tx  = mean(tx, na.rm = TRUE),
                                 mean_tn  = mean(tn, na.rm = TRUE),
                                 sum_prcp = sum(prcp, na.rm = TRUE)) %>%
                dplyr::mutate(proporcion_faltantes_tmax = 1 - cantidad_tmax/cantidad_datos,
                              proporcion_faltantes_tmin = 1 - cantidad_tmin/cantidad_datos,
                              proporcion_faltantes_prcp = 1 - cantidad_prcp/cantidad_datos) %>%
                dplyr::mutate(mean_tx  = dplyr::if_else(proporcion_faltantes_tmax > umbral.faltantes, as.double(NA), mean_tx),
                              mean_tn  = dplyr::if_else(proporcion_faltantes_tmin > umbral.faltantes, as.double(NA), mean_tn),
                              sum_prcp = dplyr::if_else(proporcion_faltantes_prcp > umbral.faltantes, as.double(NA), sum_prcp)) %>%
                dplyr::ungroup() %>%
                dplyr::select(station, year, season, sum_prcp, mean_tx, mean_tn) %>%
                base::as.data.frame()

            # Imputar datos faltantes (solo si es necesario)
            if (anyNA(estadisticas$sum_prcp) || anyNA(estadisticas$mean_tx) || anyNA(estadisticas$mean_tn)) {
                estadisticas.missmda   <- missMDA::imputePCA(X = dplyr::select(estadisticas, sum_prcp, mean_tx, mean_tn))
                estadisticas.imputadas <- cbind(dplyr::select(estadisticas, station, year, season), estadisticas.missmda$completeObs)
                return (estadisticas.imputadas)
            }
            return (estadisticas) # en caso que no sea necesario imputar nada, estadisticas tiene los resultados finales
        }
    )
}

create_seasonal_covariates <- function(seasonal_totals) {
    ## create seasonal total vector
    season_len = c(sum(c(31, 28, 31)), sum(c(30, 31, 30)), sum(c(31, 31, 30)), sum(c(31, 30, 31)))
    season_len_leap = c(sum(c(31, 29, 31)), sum(c(30, 31, 30)), sum(c(31, 31, 30)), sum(c(31, 30, 31)))

    season_rainfall_covariates <- list(c(), c(), c(), c())
    season_max_temp_covariates <- list(c(), c(), c(), c())
    season_min_temp_covariates <- list(c(), c(), c(), c())

    # ST1 = ST2 = ST3 = ST4 = c() SMX1 = SMX2 = SMX3 = SMX4 = c() SMN1 = SMN2 = SMN3 = SMN4 = c()

    for (yr in unique(seasonal_totals$year)) {
        for (season in unique(seasonal_totals$season)) {
            season_data <- seasonal_totals[seasonal_totals$year == yr & seasonal_totals$season == season, ]
            season_len_data <- if (lubridate::leap_year(yr)) season_len_leap else season_len

            values_before <- if (season > 1)
                rep(0, times = sum(season_len_data[1:(season - 1)])) else NULL
            values_after <- if (season < 4)
                rep(0, times = sum(season_len_data[(season + 1):4])) else NULL

            season_rainfall_covariates[[season]] <- c(season_rainfall_covariates[[season]], values_before, rep(season_data$sum_prcp, times = season_len_data[season]), values_after)

            season_max_temp_covariates[[season]] <- c(season_max_temp_covariates[[season]], values_before, rep(season_data$mean_tx, times = season_len_data[season]), values_after)

            season_min_temp_covariates[[season]] <- c(season_min_temp_covariates[[season]], values_before, rep(season_data$mean_tn, times = season_len_data[season]), values_after)
        }
    }

    return(list(tx = season_max_temp_covariates, tn = season_min_temp_covariates, prcp = season_rainfall_covariates))
}


# rm(season, values_after, values_before, yr, season_data, season_len_data, d_montot, d_avg_montot)
