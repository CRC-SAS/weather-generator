
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
                dplyr::mutate(season = lubridate::quarter(date, fiscal_start = 12), year = lubridate::year(date))

            estadisticas <- datos_stn_id %>%
                dplyr::group_by(station_id, year, season) %>%
                dplyr::summarise(cantidad_datos = dplyr::n(),
                                 cantidad_tmax = sum(ifelse(is.na(tmax), 0, 1)),
                                 cantidad_tmin = sum(ifelse(is.na(tmin), 0, 1)),
                                 cantidad_prcp = sum(ifelse(is.na(prcp), 0, 1)),
                                 mean_tmax = mean(tmax, na.rm = TRUE),
                                 mean_tmin = mean(tmin, na.rm = TRUE),
                                 sum_prcp  = sum(prcp, na.rm = TRUE)) %>%
                dplyr::mutate(proporcion_faltantes_tmax = 1 - cantidad_tmax/cantidad_datos,
                              proporcion_faltantes_tmin = 1 - cantidad_tmin/cantidad_datos,
                              proporcion_faltantes_prcp = 1 - cantidad_prcp/cantidad_datos) %>%
                dplyr::mutate(mean_tmax = dplyr::if_else(proporcion_faltantes_tmax > umbral_faltantes, as.double(NA), mean_tmax),
                              mean_tmin = dplyr::if_else(proporcion_faltantes_tmin > umbral_faltantes, as.double(NA), mean_tmin),
                              sum_prcp  = dplyr::if_else(proporcion_faltantes_prcp > umbral_faltantes, as.double(NA), sum_prcp)) %>%
                dplyr::ungroup() %>%
                dplyr::select(station_id, year, season, sum_prcp, mean_tmax, mean_tmin)

            # Imputar datos faltantes (solo si es necesario)
            if (anyNA(estadisticas$sum_prcp) || anyNA(estadisticas$mean_tmax) || anyNA(estadisticas$mean_tmin)) {
                estadisticas_as_df     <- estadisticas %>% base::as.data.frame()
                estadisticas_missmda   <- missMDA::imputePCA(X = dplyr::select(estadisticas_as_df, sum_prcp, mean_tmax, mean_tmin))
                estadisticas_imputadas <- cbind(dplyr::select(estadisticas, station_id, year, season), estadisticas_missmda$completeObs)
                return (estadisticas_imputadas %>% tibble::as_tibble())
            }
            return (estadisticas) # en caso que no sea necesario imputar nada, estadisticas tiene los resultados finales
        }
    )
}
