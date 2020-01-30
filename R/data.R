

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
climate.tibble <- function(date, station_id, tmax, tmin, prcp, ...) {
    return (tibble::tibble(date = as.Date(date),
                           station_id = as.integer(station_id),
                           tmax = as.numeric(tmax),
                           tmin = as.numeric(tmin),
                           prcp = as.numeric(prcp), ...))
}

#' @title Transform object to climate tibble
#' @description Transform object to climate tibble
#' @export
as.climate.tibble <- function(object, map_cols = list(date = "date", station_id = "station_id",
                                                      tmax = "tmax", tmin = "tmin", prcp = "prcp")) {
    req_cols <- list(date = "date", station_id = "station_id",
                     tmax = "tmax", tmin = "tmin", prcp = "prcp")
    lst_cols <- setdiff(names(req_cols), names(map_cols))
    map_cols <- c(map_cols, req_cols[lst_cols])

    climate <- object %>% tibble::as_tibble() %>%
        dplyr::rename(date = !!rlang::sym(map_cols$date),
                      station_id = !!rlang::sym(map_cols$station_id),
                      tmax = !!rlang::sym(map_cols$tmax),
                      tmin = !!rlang::sym(map_cols$tmin),
                      prcp = !!rlang::sym(map_cols$prcp)) %>%
        dplyr::mutate(date = as.Date(date),
                      station_id = as.integer(station_id),
                      tmax = as.numeric(tmax),
                      tmin = as.numeric(tmin),
                      prcp = as.numeric(prcp)) %>%
        tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days")) %>%
        exclude_incomplete_years() # solo se pueden tomar anhos completos!!
    return (climate)
}


#' @title Create stations simple feature (sf)
#' @description Create stations simple feature (sf).
#' @export
stations.sf <- function(station_id, longitude, latitude, crs, ...) {
    stations <- tibble::tibble(station_id = as.integer(station_id),
                               longitude = as.numeric(longitude),
                               latitude = as.numeric(latitude), ...) %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = crs)
    return (stations)
}

#' @title Transform object to stations simple feature (sf)
#' @description Transform object to stations simple feature (sf).
#' @export
as.stations.sf <- function(object, crs, map_cols = list(station_id = "station_id",
                                                        longitude = "lon_dec", latitude = "lat_dec")) {
    req_cols <- list(station_id = "station_id", longitude = "lon_dec", latitude = "lat_dec")
    lst_cols <- setdiff(names(req_cols), names(map_cols))
    map_cols <- c(map_cols, req_cols[lst_cols])

    stations <- object %>% tibble::as_tibble() %>%
        dplyr::rename(station_id = !!rlang::sym(map_cols$station_id),
                      longitude = !!rlang::sym(map_cols$longitude),
                      latitude = !!rlang::sym(map_cols$latitude) ) %>%
        dplyr::mutate(station_id = as.integer(station_id),
                      longitude = as.numeric(longitude),
                      latitude = as.numeric(latitude)) %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = crs)
    return (stations)
}


#' @title Create simulation locations simple feature (sf)
#' @description Create simulation locations simple feature (sf).
#' @export
sim_loc.sf <- function(station_id, longitude, latitude, crs, ...) {
    stations <- tibble::tibble(station_id = as.integer(station_id),
                               longitude = as.numeric(longitude),
                               latitude = as.numeric(latitude), ...) %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = crs)
    return (stations)
}

#' @title Transform object to simulation locations simple feature (sf)
#' @description Transform object to simulation locations simple feature (sf).
#' @export
as.sim_loc.sf <- function(object, crs, map_cols = list(station_id = "station_id",
                                                       longitude = "longitude", latitude = "latitude")) {
    req_cols <- list(station_id = "station_id", longitude = "longitude", latitude = "latitude")
    lst_cols <- setdiff(names(req_cols), names(map_cols))
    map_cols <- c(map_cols, req_cols[lst_cols])

    stations <- object %>% tibble::as_tibble() %>%
        dplyr::rename(station_id = !!rlang::sym(map_cols$station_id),
                      longitude = !!rlang::sym(map_cols$longitude),
                      latitude = !!rlang::sym(map_cols$latitude) ) %>%
        dplyr::mutate(station_id = as.integer(station_id),
                      longitude = as.numeric(longitude),
                      latitude = as.numeric(latitude)) %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = crs)
    return (stations)
}


check.object <- function(object, obj.name, obj.class, obj.cols) {

    if (!all(class(object) %in% obj.class)) {
        warning (glue::glue("{obj.name} must be a {if ('sf' %in% obj.class) 'sf'} tibble"), call. = FALSE, immediate. = TRUE)
        stop ("Entry data aren't in the correct format!", call. = FALSE)
    }

    if (!all(names(obj.cols) %in% names(object))) {
        warning (paste(glue::glue("{obj.name} must have the following columns:"),
                       glue::glue("{toString(lapply(names(obj.cols), function(c) glue::glue('{c} ({obj.cols[[c]]})')))}")),
                 call. = FALSE, immediate. = TRUE)
        stop ("Entry data aren't in the correct format!", call. = FALSE)
    }

    for (c in names(obj.cols)) {
        if (class(dplyr::pull(object,!!c))[[1]] != obj.cols[[c]]) {
            warning (glue::glue("{obj.name} column named {c} must be {obj.cols[[c]]}"), call. = FALSE, immediate. = TRUE)
            stop ("Entry data aren't in the correct format!", call. = FALSE)
        }
    }

}

check.ends.with.columns <- function(object, obj.name, obj.cols) {

    for (c in names(obj.cols)) {
        if (sum(endsWith(names(object), c)) != 1) {
            warning (glue::glue("{obj.name} must have exactly one column ending with: {c}!"),
                     call. = FALSE, immediate. = TRUE)
            stop ("Entry data aren't in the correct format!", call. = FALSE)
        }
    }

    for (c in names(obj.cols)) {
        if (class(dplyr::pull(dplyr::select(object, dplyr::ends_with(c))))[[1]] != obj.cols[[c]]) {
            warning (glue::glue("{obj.name} column named {c} must be {obj.cols[[c]]}"), call. = FALSE, immediate. = TRUE)
            stop ("Entry data aren't in the correct format!", call. = FALSE)
        }
    }

}

check.fit.input.data <- function(climate, stations, seasonal.climate) {

    climate.columns <- c(date = "Date", station_id = "integer", tmax = "numeric", tmin = "numeric", prcp = "numeric")
    stations.columns <- c(station_id = "integer", geometry = "sfc_POINT")
    seasonal.invariant.columns <- c(station_id = "integer", year = "numeric", season = "numeric")
    seasonal.ends.with.columns <- c(tmax = "numeric", tmin = "numeric", prcp = "numeric")

    check.object(climate, "climate", c("tbl_df", "tbl", "data.frame"), climate.columns)
    check.object(stations, "stations", c("sf", "tbl_df", "tbl", "data.frame"), stations.columns)

    if (!is.null(seasonal.climate)) {
        check.object(seasonal.climate, "seasonal.climate", c("tbl_df", "tbl", "data.frame"), seasonal.invariant.columns)
        check.ends.with.columns(seasonal.climate, "seasonal.climate", seasonal.ends.with.columns)
    }

}

check.simulation.input.data <- function(simulation.locations, seasonal.climate) {

    sim.loc.columns <- c(geometry = "sfc_POINT")
    seasonal.invariant.columns <- c(station_id = "integer", year = "numeric", season = "numeric")
    seasonal.ends.with.columns <- c(tmax = "numeric", tmin = "numeric", prcp = "numeric")

    check.object(simulation.locations, "simulation_locations", c("sf", "tbl_df", "tbl", "data.frame"), sim.loc.columns)

    if (!is.null(seasonal.climate)) {
        check.object(seasonal.climate, "seasonal_climate", c("tbl_df", "tbl", "data.frame"), seasonal.invariant.columns)
        check.ends.with.columns(seasonal.climate, "seasonal.climate", seasonal.ends.with.columns)
    }

}

check.points.to.extract <- function(points_to_extract) {

    points.columns <- c(point_id = "integer", geometry = "sfc_POINT")
    check.object(points_to_extract, "points_to_extract", c("sf", "tbl_df", "tbl", "data.frame"), points.columns)

}


#' @title Transform netcdf4 file to tibble
#' @description Transform netcdf4 file to tibble.
#' @export
netcdf.as.tibble <- function(netcdf_filename, points_to_extract, points_id_column = "point_id") {

    # Determinar variables y cantidad de realizaciones
    netcdf_file          <- ncdf4::nc_open(filename = netcdf_filename)
    variables            <- names(netcdf_file$var)
    numero_realizaciones <- netcdf_file$dim$realization$len
    coord_ref_system     <- ncdf4::ncatt_get(netcdf_file,0,"CRS")$value
    ncdf4::nc_close(netcdf_file)
    rm(netcdf_file)

    # Transform points to the correct crs
    points <- points_to_extract %>%
        sf::st_transform(crs = coord_ref_system) %>%
        dplyr::select(point_id = !!points_id_column)

    # Verificar columnas del objeto points_to_extract
    glmwgen:::check.points.to.extract(points)

    # Determinar posiciones de cada estacion (fila/columna)
    first_brick <- raster::brick(netcdf_filename, varname = variables[1], lvar = 4, level = 1,
                                 stopIfNotEqualSpaced=FALSE) %>% dplyr::first()
    cell_of_pnt <- raster::extract(x = first_brick, y = points, cellnumbers = TRUE, df = TRUE) %>%
        dplyr::select(cells) %>% dplyr::rename(cell = cells) %>%
        dplyr::bind_cols(points) %>% tidyr::drop_na(cell)

    # Obtener los datos de todas las variables y realizaciones en un data frame
    datos_simulaciones <- purrr::pmap_dfr(
        .l = purrr::cross2(variables, seq_len(numero_realizaciones)) %>% purrr::transpose(),
        .f = function(variable, numero_realizacion) {
            brick_variable <- raster::brick(netcdf_filename, varname = variable, lvar = 4, level = numero_realizacion)
            points_data  <- raster::extract(x = brick_variable, y = dplyr::pull(cell_of_pnt, cell), df = TRUE) %>%
                dplyr::mutate(variable = variable, realizacion = numero_realizacion, point_id = dplyr::pull(cell_of_pnt, point_id))
        }
    ) %>% tidyr::pivot_longer(cols = starts_with("X"), names_to = "fecha_string", values_to = "valor") %>%
        dplyr::mutate(fecha = as.Date(fecha_string, format ="X%Y.%m.%d")) %>%
        dplyr::select(point_id, realizacion, fecha, variable, valor) %>%
        tidyr::pivot_wider(names_from = "variable", values_from = "valor") %>%
        dplyr::arrange(point_id, realizacion, fecha) %>%
        dplyr::mutate(tipo_dia = factor(ifelse(as.logical(prcp), 'Lluvioso', 'Seco'), levels = c('Lluvioso', 'Seco')))

    return (datos_simulaciones)

}

