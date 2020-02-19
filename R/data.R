
# This function returns True when the seasonal covariates aren't related
# to any particular station or simulation point, and so, when the seasonal
# covariates must be repeated for any simulation point or station instead
# of being interpolated!!!
repeat_seasonal_covariates <- function(seasonal_covariates) {

    # If seasonal_covariates have a column named station_id,
    # it don't must be repeated for each simulation point
    if("station_id" %in% colnames(seasonal_covariates))
        return (FALSE)

    return (TRUE)

}


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
        warning (glue::glue("{obj.name} must be a {if ('sf' %in% obj.class) 'sf tibble' else 'tibble'}"), call. = FALSE, immediate. = TRUE)
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

check.one.of.columns <- function(object, obj.name, obj.cols) {

    if (!any(names(obj.cols) %in% names(object))) {
        warning (glue::glue("{obj.name} must have one of this columns {paste(names(obj.cols), collapse = ', ')}!"),
                 call. = FALSE, immediate. = TRUE)
        stop ("Entry data aren't in the correct format!", call. = FALSE)
    }

    for (c in names(obj.cols)) {
        if (c %in% names(object)) {
            if (class(dplyr::pull(object,!!c))[[1]] != obj.cols[[c]]) {
                warning (glue::glue("{obj.name} column named {c} must be {obj.cols[[c]]}"), call. = FALSE, immediate. = TRUE)
                stop ("Entry data aren't in the correct format!", call. = FALSE)
            }
        }
    }

}

check.fit.input.data <- function(climate, stations, seasonal.covariates) {

    climate.columns <- c(date = "Date", station_id = "integer", tmax = "numeric", tmin = "numeric", prcp = "numeric")
    stations.columns <- c(station_id = "integer", geometry = "sfc_POINT")
    seasonal.invariant.columns <- c(station_id = "integer", year = "integer", season = "integer")
    seasonal.ends.with.columns <- c(tmax = "numeric", tmin = "numeric", prcp = "numeric")

    check.object(climate, "climate", c("tbl_df", "tbl", "data.frame"), climate.columns)
    check.object(stations, "stations", c("sf", "tbl_df", "tbl", "data.frame"), stations.columns)

    if (!is.null(seasonal.covariates)) {
        check.object(seasonal.covariates, "seasonal_covariates", c("tbl_df", "tbl", "data.frame"), seasonal.invariant.columns)
        check.ends.with.columns(seasonal.covariates, "seasonal_covariates", seasonal.ends.with.columns)
    }

}

check.simulation.input.data <- function(simulation.locations, seasonal.covariates) {

    repeatable_seasonal_cov <- gamwgen:::repeat_seasonal_covariates(seasonal.covariates)

    sim.loc.columns <- c(geometry = "sfc_POINT")
    if (repeatable_seasonal_cov)
        seasonal.invariant.columns <- c(year = "integer", season = "integer")
    if (!repeatable_seasonal_cov)
        seasonal.invariant.columns <- c(station_id = "integer", year = "integer", season = "integer")
    seasonal.ends.with.columns <- c(tmax = "numeric", tmin = "numeric", prcp = "numeric")

    check.object(simulation.locations, "simulation_locations", c("sf", "tbl_df", "tbl", "data.frame"), sim.loc.columns)

    if (!is.null(seasonal.covariates)) {
        check.object(seasonal.covariates, "seasonal_covariates", c("tbl_df", "tbl", "data.frame"), seasonal.invariant.columns)
        check.ends.with.columns(seasonal.covariates, "seasonal_covariates", seasonal.ends.with.columns)
    }

}

check.points.to.extract <- function(points_to_extract) {

    points.columns <- c(geometry = "sfc_POINT")
    check.object(points_to_extract, "points_to_extract", c("sf", "tbl_df", "tbl", "data.frame"), points.columns)
    points.one.of.columns <- c(point_id = "integer", station_id = "integer")
    check.one.of.columns(points_to_extract, "points_to_extract", points.one.of.columns)

}


#' @title Extract specific points from netcdf4, as tibble
#' @description Extract specific points from a netcdf4 file, as a tibble object.
#' @import dplyr
#' @import magrittr
#' @export
netcdf.extract.points.as.tibble <- function(netcdf_filename, points_to_extract) {

    # Se carga el paquete magrittr para poder usar %T>% y apagar el warning lanzado por dplyr::one_of
    suppressPackageStartupMessages(library("magrittr"))
    old_warn_value <- getOption("warn")

    # Verificar columnas del objeto points_to_extract
    gamwgen:::check.points.to.extract(points_to_extract)

    # Determinar variables y cantidad de realizaciones
    netcdf_file          <- ncdf4::nc_open(filename = netcdf_filename)
    variables            <- names(netcdf_file$var)
    numero_realizaciones <- netcdf_file$dim$realization$len
    coord_ref_system     <- ncdf4::ncatt_get(netcdf_file,0,"CRS")$value
    ncdf4::nc_close(netcdf_file)
    rm(netcdf_file)

    # Transform points to the correct crs
    points <- points_to_extract %>%
        sf::st_transform(crs = coord_ref_system) %T>%
        { base::options(warn=-1) } %>%
        dplyr::select(dplyr::one_of("station_id", "point_id")) %T>%
        { base::options(warn=old_warn_value) } %>%
        dplyr::mutate(ID = 1:dplyr::n())

    # Determinar posiciones de cada estacion (fila/columna)
    first_brick <- raster::brick(netcdf_filename, varname = variables[1], lvar = 4, level = 1,
                                 stopIfNotEqualSpaced=FALSE) %>% dplyr::first()
    cell_of_pnt <- raster::extract(x = first_brick, y = points, cellnumbers = TRUE, df = TRUE) %>%
        dplyr::select(cells) %>% dplyr::rename(cell = cells) %>%
        dplyr::bind_cols(points) %>% tidyr::drop_na(cell)

    # Obtener los datos de todas las variables y realizaciones en un data frame
    datos_simulaciones_wide <- purrr::pmap_dfr(
        .l = purrr::cross2(variables, seq_len(numero_realizaciones)) %>% purrr::transpose(),
        .f = function(variable, numero_realizacion) {
            brick_variable <- raster::brick(netcdf_filename, varname = variable, lvar = 4, level = numero_realizacion)
            points_data  <- raster::extract(x = brick_variable, y = dplyr::pull(cell_of_pnt, cell), df = TRUE) %>%
                dplyr::mutate(variable = variable, realization = numero_realizacion, ctrl_id = dplyr::pull(cell_of_pnt, ID))
        })

    # Reestructurar datos_simulaciones
    datos_simulaciones_large <- datos_simulaciones_wide %>%
        tidyr::pivot_longer(cols = starts_with("X"), names_to = "date_string", values_to = "valor") %>%
        dplyr::mutate(date = as.Date(date_string, format ="X%Y.%m.%d")) %>% dplyr::select(-date_string) %>%
        tidyr::pivot_wider(names_from = "variable", values_from = "valor")

    # Modificaciones finales (se agrega los ids recibidos en points_to_extract)
    datos_simulaciones <- datos_simulaciones_large %>%
        dplyr::inner_join(points, by = "ID") %>% dplyr::select(-geometry) %>%
        dplyr::mutate(type_day = base::factor(x = ifelse(as.logical(prcp), 'Wet', 'Dry'),
                                              levels = c('Wet', 'Dry'))) %T>%
        { base::options(warn=-1) } %>%
        dplyr::select(dplyr::one_of("station_id", "point_id"), realization, date,
                      tmax, tmin, prcp, type_day) %>%
        dplyr::arrange_at(dplyr::vars(realization, dplyr::one_of("station_id", "point_id"), date)) %T>%
        { base::options(warn=old_warn_value) }

    return (datos_simulaciones)

}


#' @title Transform netcdf4 file to tibble object
#' @description Transform netcdf4 file to tibble object.
#' @export
netcdf.as.tibble <- function(netcdf_filename, na.rm = T) {

    result <- tidync::tidync(netcdf_filename) %>%
        tidync::hyper_tibble(na.rm = na.rm) %>%
        dplyr::mutate(time = lubridate::as_date(time),
                      type_day = base::factor(x = ifelse(as.logical(prcp), 'Wet', 'Dry'),
                                              levels = c('Wet', 'Dry'))) %>%
        dplyr::rename(date = time,
                      longitude = projection_x_coordinate,
                      latitude = projection_y_coordinate) %>%
        dplyr::select(realization, date, tmax, tmin, prcp, type_day, longitude, latitude)

    return (result)

}


#' @title Transform netcdf4 file to sf object
#' @description Transform netcdf4 file to sf object.
#' @export
netcdf.as.sf <- function(netcdf_filename, na.rm = T) {

    nc_proj4string <- ncmeta::nc_att(netcdf_filename, "NC_GLOBAL", "CRS")$value$CRS
    nc_crs  <- sf::st_crs(nc_proj4string)

    result <- tidync::tidync(netcdf_filename) %>%
        tidync::hyper_tibble(na.rm = na.rm) %>%
        dplyr::mutate(time = lubridate::as_date(time),
                      type_day = base::factor(x = ifelse(as.logical(prcp), 'Wet', 'Dry'),
                                              levels = c('Wet', 'Dry')),
                      longitude = projection_x_coordinate,
                      latitude = projection_y_coordinate) %>%
        dplyr::rename(date = time) %>%
        sf::st_as_sf(coords = c('projection_x_coordinate', 'projection_y_coordinate'),
                     crs = nc_crs) %>%
        dplyr::select(realization, date, tmax, tmin, prcp, type_day, longitude, latitude)

    return (result)

}


#' @title Extract specific points from netcdf4, as sf
#' @description Extract specific points from a netcdf4 file, as a sf object.
#' @import dplyr
#' @import magrittr
#' @export
netcdf.extract.points.as.sf <- function(netcdf_filename, points_to_extract) {

    # Se carga el paquete magrittr para poder usar %T>% y apagar el warning lanzado por dplyr::one_of
    suppressPackageStartupMessages(library("magrittr"))
    old_warn_value <- getOption("warn")

    # Verificar columnas del objeto points_to_extract
    gamwgen:::check.points.to.extract(points_to_extract)

    # Determinar la proyección de las coordenadas en el NetCDF
    nc_proj4string <- ncmeta::nc_att(netcdf_filename, "NC_GLOBAL", "CRS")$value$CRS
    nc_crs  <- sf::st_crs(nc_proj4string)

    # Extraer datos del NetCDF (todos los datos)
    netcdf_as_sf <- gamwgen::netcdf.as.sf(netcdf_filename)

    # Extraer las coordenadas de las ubiciones en el netcdf
    netcdf_points <- netcdf_as_sf %>%
        dplyr::filter(realization == 1) %>%
        sf::st_drop_geometry() %>%
        dplyr::distinct(longitude, latitude) %>%
        dplyr::mutate(lon = longitude, lat = latitude) %>%
        sf::st_as_sf(coords = c('lon', 'lat'), crs = nc_crs) %>%
        dplyr::mutate(ID = 1:dplyr::n())

    # Obtener la menor distancia entre las ubicaciones extraídas
    min_distance <- 10000
    if (nrow(netcdf_points) > 1)
        min_distance <- purrr::map_dfr(
            .x = netcdf_points %>% dplyr::pull(ID),
            .f = function(id, all_pts) {
                target_point  <- all_pts %>% dplyr::filter(ID == id)
                remaining_pts <- all_pts %>% dplyr::filter(ID != id)
                remianing_nid <- sf::st_nearest_feature(target_point, remaining_pts)
                nearest_point <- remaining_pts %>% dplyr::slice(remianing_nid)
                dist_nearest  <- sf::st_distance(target_point, nearest_point) %>% as.integer()
                return (tibble::tibble(ID = id, nID = nid, distance = dist_nearest))
            }, all_pts = netcdf_points) %>%
            dplyr::pull(distance) %>% min()

    # Transformar cada punto extraído del netcdf en un poligono cuadrado de lado igual a
    # min_distance/2, esto para poder intersectar estos poligonos con los puntos a extraer!!
    netcdf_as_sf_polygonized <- netcdf_as_sf %>%
        sf::st_buffer(dist = min_distance/2, endCapStyle="SQUARE")

    # Transform points to the correct crs
    points_to_extract <- points_to_extract %>%
        sf::st_transform(crs = nc_crs) %T>%
        { base::options(warn=-1) } %>%
        dplyr::select(-dplyr::one_of("longitude", "latitude")) %T>%
        { base::options(warn=old_warn_value) }

    # Se seleccionan solo los resultados deseados
    resulting_data <- netcdf_as_sf_polygonized %>%
        sf::st_join(points_to_extract, left = FALSE) %>% sf::st_drop_geometry() %>% tibble::as_tibble() %T>%
        { base::options(warn=-1) } %>%
        dplyr::select(dplyr::one_of("station_id", "point_id"), realization, date, tmax, tmin, prcp, type_day,
                      longitude, latitude) %>%
        dplyr::arrange_at(dplyr::vars(realization, dplyr::one_of("station_id", "point_id"), date)) %T>%
        { base::options(warn=old_warn_value) } %>%
        dplyr::mutate(projection_x_coordinate = longitude, projection_y_coordinate = latitude) %>%
        sf::st_as_sf(coords = c("projection_x_coordinate", "projection_y_coordinate"), crs = nc_crs)

    return (resulting_data)

}

