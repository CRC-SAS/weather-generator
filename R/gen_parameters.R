
start_climatology_month_day <- function(model, simulation_points, month, day) {

    start_climatology <- model$start_climatology %>%
        dplyr::filter(month == lubridate::month(start_date-1),
                      day == lubridate::day(start_date-1)) %>%
        dplyr::mutate(prcp_occ = as.integer(prcp > 0))

    start_climatology_sf <- simulation_points %>%
        dplyr::left_join(start_climatology, by = "station_id") %>%
        dplyr::select(prcp_occ, tmax, tmin)

    return (start_climatology_sf)

}

generate_random_temperatura_noise <- function(residuals_statistics, simulation_dates) {
    ruidos_aleatorios <- residuals_statistics %>%
        dplyr::filter(month %in% unique(simulation_dates$month)) %>%
        purrr::pmap_dfr(
            function(station_id, month, tipo_dia, sd.tmax_residuals, sd.tmin_residuals,
                     mean.tmax_residuals, mean.tmin_residuals, cov.residuals,
                     var.tmax_residuals, var.tmin_residuals) {
                ruidos <- MASS::mvrnorm(n = nrow(simulation_dates %>% dplyr::filter(month == !!month)),
                                        mu = c(mean.tmax_residuals, mean.tmin_residuals),
                                        Sigma = matrix(c(var.tmax_residuals, cov.residuals, cov.residuals, var.tmin_residuals), 2, 2),
                                        empirical = TRUE) %>%
                    magrittr::set_colnames(c("tmax_noise", "tmin_noise"))
                tibble::tibble(station_id = as.integer(station_id),
                               date       = simulation_dates %>% dplyr::filter(month == !!month) %>% dplyr::pull(date),
                               tipo_dia   = tipo_dia,
                               tmax_noise = ruidos[,"tmax_noise"],
                               tmin_noise = ruidos[,"tmin_noise"])
            })
}

generate_residuals_statistics <- function(models_residuals) {

    residuals_statistics <- models_residuals %>% tidyr::drop_na() %>%
        dplyr::select(station_id, date, tmax_residuals, tmin_residuals, tipo_dia) %>%
        dplyr::mutate(month = lubridate::month(date)) %>%
        dplyr::group_by(station_id, month) %>%
        dplyr::mutate(sd.tmax_residuals = stats::sd(tmax_residuals, na.rm = T),
                      sd.tmin_residuals = stats::sd(tmin_residuals, na.rm = T)) %>%
        dplyr::group_by(tipo_dia, add = T) %>%
        dplyr::mutate(mean.tmax_residuals = mean(tmax_residuals, na.rm = T),
                      mean.tmin_residuals = mean(tmin_residuals, na.rm = T),
                      cov.residuals       = stats::cov(tmax_residuals, tmin_residuals, use = "pairwise.complete.obs"),
                      var.tmax_residuals  = stats::var(tmax_residuals, na.rm = T),
                      var.tmin_residuals  = stats::var(tmin_residuals, na.rm = T)) %>%
        dplyr::select(-date, -tmax_residuals, -tmin_residuals) %>%
        dplyr::ungroup() %>% dplyr::distinct()

    return (residuals_statistics)
}

generate_month_params <- function(residuals, observed_climate, stations) {
    month_params <- purrr::map(
        .x = 1:12,
        .f = function(m) {

            ## OBS:
            # Estamos trabajando con coordenadas métricas, sin embargo,
            # el paquete RandomFields trabaja mejor en metros, por la tanto,
            # pasamos distance_matrix a kilometros, dividiendo por 1000!!
            distance_matrix <- glmwgen:::make_distance_matrix(stations) / 1000

            month_residuals <- residuals %>%
                dplyr::filter(lubridate::month(date) == m) %>%
                dplyr::mutate(., year = lubridate::year(date)) %>%
                dplyr::arrange(date, station_id)

            month_climate   <- observed_climate %>%
                dplyr::filter(lubridate::month(date) == m) %>%
                dplyr::mutate(., year = lubridate::year(date)) %>%
                dplyr::select(date, year, station_id, tmax, tmin, prcp, tipo_dia) %>%
                tidyr::drop_na(.)

            # Crear matriz de residuos de temperatura maxima en formato ancho
            tmax_residuals_matrix_dry_days <- month_residuals %>%
                dplyr::filter(., tipo_dia == 'Seco') %>%
                dplyr::select(., station_id, date, tmax_residuals) %>%
                tidyr::spread(., key = station_id, value = tmax_residuals) %>%
                dplyr::select(., colnames(distance_matrix))
            tmax_residuals_matrix_wet_days <- month_residuals %>%
                dplyr::filter(., tipo_dia == 'Lluvioso') %>%
                dplyr::select(., station_id, date, tmax_residuals) %>%
                tidyr::spread(., key = station_id, value = tmax_residuals) %>%
                dplyr::select(., colnames(distance_matrix))

            # Crear matriz de residuos de temperatura minima en formato ancho
            tmin_residuals_matrix_dry_days <- month_residuals %>%
                dplyr::filter(., tipo_dia == 'Seco') %>%
                dplyr::select(., station_id, date, tmin_residuals) %>%
                tidyr::spread(., key = station_id, value = tmin_residuals) %>%
                dplyr::select(., colnames(distance_matrix))
            tmin_residuals_matrix_wet_days <- month_residuals %>%
                dplyr::filter(., tipo_dia == 'Lluvioso') %>%
                dplyr::select(., station_id, date, tmin_residuals) %>%
                tidyr::spread(., key = station_id, value = tmin_residuals) %>%
                dplyr::select(., colnames(distance_matrix))

            # Crear matriz de residuos de ocurrencia de lluvia en formato ancho
            probit_residuals_matrix <- month_residuals %>%
                dplyr::select(., station_id, date, prcp_occ_residuals) %>%
                tidyr::spread(., key = station_id, value = prcp_occ_residuals) %>%
                dplyr::select(., colnames(distance_matrix))

            # Matrix de datos observados de temperatura maxima en formato ancho
            tmax_matrix_dry_days <- month_climate %>%
                dplyr::filter(., tipo_dia == 'Seco') %>%
                dplyr::select(., station_id, date, tmax) %>%
                tidyr::spread(., key = station_id, value = tmax) %>%
                dplyr::select(., colnames(distance_matrix))
            tmax_matrix_wet_days <- month_climate %>%
                dplyr::filter(., tipo_dia == 'Lluvioso') %>%
                dplyr::select(., station_id, date, tmax) %>%
                tidyr::spread(., key = station_id, value = tmax) %>%
                dplyr::select(., colnames(distance_matrix))

            # Matrix de datos observados de temperatura minima en formato ancho
            tmin_matrix_dry_days <- month_climate %>%
                dplyr::filter(., tipo_dia == 'Seco') %>%
                dplyr::select(., station_id, date, tmin) %>%
                tidyr::spread(., key = station_id, value = tmin) %>%
                dplyr::select(., colnames(distance_matrix))
            tmin_matrix_wet_days <- month_climate %>%
                dplyr::filter(., tipo_dia == 'Lluvioso') %>%
                dplyr::select(., station_id, date, tmin) %>%
                tidyr::spread(., key = station_id, value = tmin) %>%
                dplyr::select(., colnames(distance_matrix))

            n_stations <- length(unique(stations))

            # Crear matrices de covarianza
            prcp_cor <- stats::cor(probit_residuals_matrix, use = "pairwise.complete")
            prcp_vario <- stats::var(probit_residuals_matrix, use = "pairwise.complete") * (1 - prcp_cor)
            tmax_vario_dry <- stats::cov(tmax_residuals_matrix_dry_days, use = "pairwise.complete.obs")
            tmax_vario_wet <- stats::cov(tmax_residuals_matrix_wet_days, use = "pairwise.complete.obs")
            tmin_vario_dry <- stats::cov(tmin_residuals_matrix_dry_days, use = "pairwise.complete.obs")
            tmin_vario_wet <- stats::cov(tmin_residuals_matrix_wet_days, use = "pairwise.complete.obs")

            # Assert that the column and row names of the variograms equal the ones of the distance matrix.
            stopifnot(all(colnames(tmax_vario_dry) == colnames(distance_matrix)), all(rownames(tmax_vario_dry) == rownames(distance_matrix)))
            stopifnot(all(colnames(tmax_vario_wet) == colnames(distance_matrix)), all(rownames(tmax_vario_wet) == rownames(distance_matrix)))
            stopifnot(all(colnames(tmin_vario_dry) == colnames(distance_matrix)), all(rownames(tmin_vario_dry) == rownames(distance_matrix)))
            stopifnot(all(colnames(tmin_vario_wet) == colnames(distance_matrix)), all(rownames(tmin_vario_wet) == rownames(distance_matrix)))

            # Estimar variograma de precipitacion
            prcp_params <- stats::optim(par = c(0.01, 1, max(distance_matrix)), f = glmwgen:::partially_apply_LS(prcp_vario, distance_matrix))$par
            prcp_params[prcp_params < 0] <- 0
            prcp_params <- c(0, 1, prcp_params[3])

            # Variograma de temperatura máxima para los días secos
            sill_initial_value_dry <- mean(var(tmax_matrix_dry_days, na.rm = T))
            tmax_params_dry <- stats::optim(par = c(sill_initial_value_dry, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmax_vario_dry, distance_matrix, base_p = c(0)))$par
            tmax_params_dry <- c(0, tmax_params_dry)
            # Variograma de temperatura máxima para los días lluviosos
            sill_initial_value_wet <- mean(var(tmax_matrix_wet_days, na.rm = T, use = 'pairwise.complete.obs'))
            tmax_params_wet <- stats::optim(par = c(sill_initial_value_wet, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmax_vario_wet, distance_matrix, base_p = c(0)))$par
            tmax_params_wet <- c(0, tmax_params_wet)

            # Variograma de temperatura mínima para los días secos
            sill_initial_value_dry <- mean(var(tmin_matrix_dry_days, na.rm = T, use = 'pairwise.complete.obs'))
            tmin_params_dry <- stats::optim(par = c(sill_initial_value_dry, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmin_vario_dry, distance_matrix, base_p = c(0)))$par
            tmin_params_dry <- c(0, tmin_params_dry)
            # Variograma de temperatura mínima para los días lluviosos
            sill_initial_value_wet <- mean(var(tmin_matrix_wet_days, na.rm = T, use = 'pairwise.complete.obs'))
            tmin_params_wet <- stats::optim(par = c(sill_initial_value_wet, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmin_vario_wet, distance_matrix, base_p = c(0)))$par
            tmin_params_wet <- c(0, tmin_params_wet)

            # Save correlation between temperature residues
            correlation <- stats::cor(month_residuals$tmax_residuals,
                                      month_residuals$tmin_residuals,
                                      use = 'complete.obs')

            # Guardar resultados
            variogram_parameters <- list(prcp = prcp_params,
                                         tmax_dry = tmax_params_dry,
                                         tmax_wet = tmax_params_wet,
                                         tmin_dry = tmin_params_dry,
                                         tmin_wet = tmin_params_wet)
            residuals_sd <- data.frame(month_residuals %>% group_by(station_id) %>%
                                           summarise(tmax = sd(tmax_residuals, na.rm = T),
                                                     tmin = sd(tmin_residuals, na.rm = T),
                                                     prcp = 1))
            rownames(residuals_sd) <- residuals_sd[, 1]
            residuals_sd <- residuals_sd[, -1]

            return(list(
                variogram_parameters = variogram_parameters,
                residuals_sd = residuals_sd,
                cov_matrix = list(
                    tmax_dry = tmax_vario_dry,
                    tmax_wet = tmax_vario_wet,
                    tmin_dry = tmin_vario_dry,
                    tmin_wet = tmin_vario_wet,
                    prcp = prcp_cor
                ),
                correlation = correlation
            ))
        }
    )
    return (month_params)
}
