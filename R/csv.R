

# Definición de función para la inserción de datos de una realización en un CSV
GuardarRealizacionEnCSV <- function(filename, numero_realizacion, tibble_with_data, avbl_cores) {

    tibble_with_data <- tibble_with_data %>%
        dplyr::select(realization = nsim,
                      tidyselect::any_of(c("station_id", "point_id")), longitude, latitude,
                      date, tmax, tmin, prcp_occ, prcp_amt, type_day)

    if (numero_realizacion == 1) {
        data.table::fwrite(tibble_with_data, file = filename, append = F, nThread = avbl_cores)
    } else {
        data.table::fwrite(tibble_with_data, file = filename, append = T, nThread = avbl_cores)
    }

}
