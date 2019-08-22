

#' @title Exclude incomplete years
#' @description Exclude incomplete years.
#' @export
exclude_incomplete_years <- function(climate_data) {
    result <- climate_data %>%
        dplyr::group_by(year = lubridate::year(date)) %>%
        dplyr::mutate(ndays = n()) %>%
        dplyr::ungroup(-year) %>%
        dplyr::filter(ndays >= 365) %>%
        dplyr::select(-year, -ndays)
    return(result)
}


#' @title Create climate tibble
#' @description Create climate tibble.
#' @export
climate.tibble <- function(date, station_id, tx, tn, prcp, ...) {
    return (tibble::tibble(date = as.Date(date), station_id = as.integer(station_id),
                           tx = as.numeric(tx), tn = as.numeric(tn),
                           prcp = as.numeric(prcp), ...))
}

#' @title Transform object to climate tibble
#' @description Transform object to climate tibble
#' @export
as.climate.tibble <- function(object, map_cols = list(date = "date", station_id = "station_id",
                                                      tx = "tx", tn = "tn", prcp = "prcp")) {
    climate <- object %>% tibble::as_tibble() %>%
        dplyr::rename(date = !!rlang::sym(map_cols$date), station_id = !!rlang::sym(map_cols$station_id),
                      tx = !!rlang::sym(map_cols$tx), tn = !!rlang::sym(map_cols$tn), prcp = !!rlang::sym(map_cols$prcp)) %>%
        dplyr::mutate(date = as.Date(date), station_id = as.integer(station_id),
                      tx = as.numeric(tx), tn = as.numeric(tn), prcp = as.numeric(prcp)) %>%
        tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days")) %>%
        exclude_incomplete_years() # solo se pueden tomar anhos completos!!
    return (climate)
}


#' @title Create stations simple feature (sf)
#' @description Create stations simple feature (sf).
#' @export
stations.sf <- function(station_id, lon_dec, lat_dec, crs, ...) {
    stations <- tibble::tibble(station_id = as.integer(station_id),
                               lon_dec = as.numeric(lon_dec), lat_dec = as.numeric(lat_dec), ...) %>%
        dplyr::mutate(latitud = lat_dec, longitud = lon_dec) %>%
        sf::st_as_sf(coords = c('longitud', 'latitud'), crs = crs)
    return (stations)
}

#' @title Transform object to stations simple feature (sf)
#' @description Transform object to stations simple feature (sf).
#' @export
as.stations.sf <- function(object, crs, map_cols = list(station_id = "station_id", lon_dec = "lon_dec", lat_dec = "lat_dec")) {
    stations <- object %>% tibble::as_tibble() %>%
        dplyr::rename(station_id = !!rlang::sym(map_cols$station_id),
                      lon_dec = !!rlang::sym(map_cols$lon_dec), lat_dec = !!rlang::sym(map_cols$lat_dec) ) %>%
        dplyr::mutate(station_id = as.integer(station_id), lon_dec = as.numeric(lon_dec), lat_dec = as.numeric(lat_dec),
                      latitud = lat_dec, longitud = lon_dec) %>%
        sf::st_as_sf(coords = c('longitud', 'latitud'), crs = crs)
    return (stations)
}


check.object <- function(object, obj.name, obj.class, obj.cols) {

    if (!all(class(object) %in% obj.class)) {
        warning (glue::glue("{obj.name} must be a {if ('sf' %in% obj.class) 'sf'} tibble"), call. = FALSE, immediate. = TRUE)
        stop ("Entry data for glmwgen::calibrate.glmwgen() aren't in the correct format!", call. = FALSE)
    }

    if (!all(names(obj.cols) %in% names(object))) {
        warning (paste(glue::glue("{obj.name} must have the following columns:"),
                       glue::glue("{toString(lapply(names(obj.cols), function(c) glue::glue('{c} ({obj.cols[[c]]})')))}")),
                 call. = FALSE, immediate. = TRUE)
        stop ("Entry data for glmwgen::calibrate.glmwgen() aren't in the correct format!", call. = FALSE)
    }

    for (c in names(obj.cols)) {
        if (class(dplyr::pull(object,!!c))[[1]] != obj.cols[[c]]) {
            warning (glue::glue("{obj.name} column named {c} must be {obj.cols[[c]]}"), call. = FALSE, immediate. = TRUE)
            stop ("Entry data for glmwgen::calibrate.glmwgen() aren't in the correct format!", call. = FALSE)
        }
    }

}

check.ends.with.columns <- function(object, obj.name, obj.cols) {

    for (c in names(obj.cols)) {
        if (sum(endsWith(names(object), c)) != 1) {
            warning (glue::glue("{obj.name} must have exactly one column ending with: {c}!"),
                     call. = FALSE, immediate. = TRUE)
            stop ("Entry data for glmwgen::calibrate.glmwgen() aren't in the correct format!", call. = FALSE)
        }
    }

    for (c in names(obj.cols)) {
        if (class(dplyr::pull(dplyr::select(object, dplyr::ends_with(c))))[[1]] != obj.cols[[c]]) {
            warning (glue::glue("{obj.name} column named {c} must be {obj.cols[[c]]}"), call. = FALSE, immediate. = TRUE)
            stop ("Entry data for glmwgen::calibrate.glmwgen() aren't in the correct format!", call. = FALSE)
        }
    }

}

check.input.data <- function(climate, stations, seasonal.climate) {

    climate.columns <- c(date = "Date", station_id = "integer", tx = "numeric", tn = "numeric", prcp = "numeric")
    stations.columns <- c(station_id = "integer", lon_dec = "numeric", lat_dec = "numeric", geometry = "sfc_POINT")
    seasonal.invariant.columns <- c(station_id = "integer", year = "numeric", season = "numeric")
    seasonal.ends.with.columns <- c(tx = "numeric", tn = "numeric", prcp = "numeric")

    check.object(climate, "climate", c("tbl_df", "tbl", "data.frame"), climate.columns)
    check.object(stations, "stations", c("sf", "tbl_df", "tbl", "data.frame"), stations.columns)

    if (!is.null(seasonal.climate)) {
        check.object(seasonal.climate, "seasonal.climate", c("tbl_df", "tbl", "data.frame"), seasonal.invariant.columns)
        check.ends.with.columns(seasonal.climate, "seasonal.climate", seasonal.ends.with.columns)
    }

}
