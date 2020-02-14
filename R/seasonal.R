
#' @title Generate seasonal covariates
#' @description Generate seasonal covariates.
#' @export
summarise_seasonal_climate <- function(datos_climaticos, umbral_faltantes = 0.2) {
    ## compute seasonal averages for tmax and tmin, and sum for prcp, with quality control and missing values imputation
    purrr::map_dfr(
        .x = unique(dplyr::pull(datos_climaticos, station_id)),
        .f = function(stn_id) {

            datos_stn_id <- datos_climaticos %>%
                dplyr::filter(station_id == stn_id) %>%
                tidyr::complete(date = base::seq.Date(from = lubridate::floor_date(min(date), "year") ,
                                                      to   = lubridate::ceiling_date(max(date), "year") - lubridate::days(1),
                                                      by   = "days"),
                                fill = list(station_id = stn_id)) %>%
                dplyr::mutate(season = lubridate::quarter(date, fiscal_start = 12),
                              year = as.integer(lubridate::year(date)), station_id = as.integer(station_id))

            estadisticas <- datos_stn_id %>%
                dplyr::group_by(station_id, year, season) %>%
                dplyr::summarise(cantidad_datos = dplyr::n(),
                                 cantidad_tmax = sum(ifelse(is.na(tmax), 0, 1)),
                                 cantidad_tmin = sum(ifelse(is.na(tmin), 0, 1)),
                                 cantidad_prcp = sum(ifelse(is.na(prcp), 0, 1)),
                                 seasonal_tmax = mean(tmax, na.rm = TRUE),
                                 seasonal_tmin = mean(tmin, na.rm = TRUE),
                                 seasonal_prcp = sum(prcp, na.rm = TRUE)) %>%
                dplyr::mutate(proporcion_faltantes_tmax = 1 - cantidad_tmax/cantidad_datos,
                              proporcion_faltantes_tmin = 1 - cantidad_tmin/cantidad_datos,
                              proporcion_faltantes_prcp = 1 - cantidad_prcp/cantidad_datos) %>%
                dplyr::mutate(seasonal_tmax = dplyr::if_else(proporcion_faltantes_tmax > umbral_faltantes, as.double(NA), seasonal_tmax),
                              seasonal_tmin = dplyr::if_else(proporcion_faltantes_tmin > umbral_faltantes, as.double(NA), seasonal_tmin),
                              seasonal_prcp = dplyr::if_else(proporcion_faltantes_prcp > umbral_faltantes, as.double(NA), seasonal_prcp)) %>%
                dplyr::ungroup() %>%
                dplyr::select(station_id, year, season, seasonal_prcp, seasonal_tmax, seasonal_tmin)

            # Imputar datos faltantes (solo si es necesario)
            if (anyNA(estadisticas$seasonal_prcp) || anyNA(estadisticas$seasonal_tmax) || anyNA(estadisticas$seasonal_tmin)) {
                estadisticas_as_df     <- estadisticas %>% base::as.data.frame()
                estadisticas_missmda   <- missMDA::imputePCA(X = dplyr::select(estadisticas_as_df, seasonal_prcp, seasonal_tmax, seasonal_tmin))
                estadisticas_imputadas <- cbind(dplyr::select(estadisticas, station_id, year, season), estadisticas_missmda$completeObs) %>%
                    dplyr::mutate(seasonal_prcp = if_else(seasonal_prcp < 0, 0, seasonal_prcp))
                return (estadisticas_imputadas %>% tibble::as_tibble())
            }
            return (estadisticas) # en caso que no sea necesario imputar nada, estadisticas tiene los resultados finales
        }
    )
}


# Definition of a function to fill gaps in the seasonal covariates data. It must be
#  used with caution given that the missing years, now filled, were not included in the fitting of the model
fill_missing_seasonal_climate <- function(seasonal_covariates) {
    ## fills missing seasonal averages for tmax and tmin, and sum for prcp, with quality control and missing values imputation
    filled.data <- purrr::map_dfr(
        .x = c('seasonal_prcp', 'seasonal_tmax', 'seasonal_tmin'),
        .f = function(variable) {

            # Selection of the variable to be filled
            seasonal.data <- seasonal_covariates %>%
                dplyr::select(station_id, year, season, variable) %>%
                tidyr::spread(., station_id, variable) %>%
                as.data.frame(.)

            # If there are missing data, the gaps will be filled.
            # Otherwise, nothing will be done to the original data
            if (anyNA(seasonal.data)) {
                # Gap filling method: missMDA
                seasonal_filled_missmda <- missMDA::imputePCA(X = seasonal.data %>%
                                                                  dplyr::select(., -c('year', 'season')))
                imputed.variable <- cbind(dplyr::select(seasonal.data, year, season), seasonal_filled_missmda$completeObs) %>%
                    tidyr::gather(station_id, value, -c('year', 'season')) %>%
                    dplyr::mutate(., variable_id = variable) %>%
                    dplyr::select(station_id, year, season, variable_id, value)

                # Prevent negative precipitation from happening
                if (variable == 'seasonal_prcp') {
                    imputed.variable %<>% dplyr::mutate(value = if_else(value < 0, 0, value))
                }

            } else {
                imputed.variable <-   seasonal_covariates %>%
                    dplyr::select(station_id, year, season, variable) %>%
                    dplyr::rename(value = variable) %>%
                    dplyr::mutate(variable_id = variable) %>%
                    dplyr::select(station_id, year, season, variable_id, value)

            }

            # return results
            return(imputed.variable)
        }
    )

    # Arrange tibble to meet the expected format
    seasonal_covariates <- filled.data %>%
        tidyr::spread(variable_id, value) %>%
        dplyr::mutate(station_id = as.integer(station_id)) %>%
        tibble::as_tibble(.)

    # Return object with filled gaps
    return(seasonal_covariates)
}

