
# Interpolación de covariables
interpolate_thresholds <- function(simulation_points, monthly_thresholds, model_stations) {

    # OBS:
    # Al interpolar las covariables, no se tiene en cuenta si el usario optó por usar ruidos
    # correlacionados espacialmente o no. Además, aunque solo sea necesario interpolar covariables
    # para una sola estación (porque para las demás ya hay covariables), de todas formas se
    # interpolan covariables para todas las estaciones, ignorando las covariables observadas
    # y utilizando en su lugar las covariables interpoladas!! Esto el funcionamiento previsto!

    # Datos estacionales para todas las localidades
    umbrales.interpolados <- monthly_thresholds %>%
        dplyr::left_join(model_stations %>% dplyr::select(station_id, longitude, latitude), by = 'station_id') %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = sf::st_crs(simulation_points))

    # Transformación a objeto sp para la interpolación
    simulation_points.sp <- simulation_points %>% sf::as_Spatial()

    # Combinaciones posibles de años y trimestres para la interpolación
    combinaciones <- purrr::transpose(
        purrr::cross2(.x = unique(monthly_thresholds$month),
                      .y = unique(monthly_thresholds$prcp_occ)))

    # Interpolación de acumulados de lluvia
    datos_interpolados <- purrr::pmap_dfr(
        .l = combinaciones,
        .f = function(month, prcp_occ) {

            interpolacion.valores.iniciales.sp <- umbrales.interpolados %>%
                dplyr::filter(month == !!month & prcp_occ == !!prcp_occ) %>%
                sf::as_Spatial()

            # Valores interpolados de ocurrencia
            datos_interpolados_max.range <-
                automap::autoKrige(max.range~1,
                                   interpolacion.valores.iniciales.sp,
                                   simulation_points.sp,
                                   debug.level = 0) %>%
                sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(grilla.simulacion)) %>%
                dplyr::select(max.range = var1.pred) %>%
                dplyr::mutate(month = !!month, prcp_occ = !!prcp_occ) %>%
                dplyr::mutate(longitude = sf::st_coordinates(geometry)[,'X'],
                              latitude = sf::st_coordinates(geometry)[,'Y']) %>%
                sf::st_drop_geometry() %>% tibble::as_tibble()


            # Valores interpolados de temperatura máxima
            datos_interpolados_min.range <-
                automap::autoKrige(min.range~1,
                                   interpolacion.valores.iniciales.sp,
                                   simulation_points.sp,
                                   debug.level = 0) %>%
                sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(grilla.simulacion)) %>%
                dplyr::select(min.range = var1.pred) %>%
                dplyr::mutate(month = !!month, prcp_occ = !!prcp_occ) %>%
                dplyr::mutate(longitude = sf::st_coordinates(geometry)[,'X'],
                              latitude = sf::st_coordinates(geometry)[,'Y']) %>%
                sf::st_drop_geometry() %>% tibble::as_tibble()


            # Crear data frame con los valores interpolados
            # Este df hay que unirlo a la grilla de simulación para cada día de la simulacion
            # para la union tenemos las coordenadas, anos y estaciones.
            datos_interpolados_range <- datos_interpolados_max.range %>%
                dplyr::left_join(datos_interpolados_min.range, by = c('month', 'prcp_occ', 'latitude', 'longitude')) %>%
                dplyr::select(month, prcp_occ, max.range, min.range, longitude, latitude)

            return(datos_interpolados_range)
        }
    )

    return (datos_interpolados)

}


# Interpolación de covariables
interpolate_covariates <- function(simulation_points, seasonal_covariates, model_stations, simulation_dates) {

    # OBS:
    # Al interpolar las covariables, no se tiene en cuenta si el usario optó por usar ruidos
    # correlacionados espacialmente o no. Además, aunque solo sea necesario interpolar covariables
    # para una sola estación (porque para las demás ya hay covariables), de todas formas se
    # interpolan covariables para todas las estaciones, ignorando las covariables observadas
    # y utilizando en su lugar las covariables interpoladas!! Esto el funcionamiento previsto!

    # Datos estacionales para todas las localidades
    datos.estacionales <- seasonal_covariates %>%
        dplyr::left_join(model_stations %>% dplyr::select(station_id, longitude, latitude), by = 'station_id') %>%
        sf::st_as_sf() %>% sf::st_transform(sf::st_crs(simulation_points))

    # Transformación a objeto sp para la interpolación
    simulation_points.sp <- simulation_points %>% sf::as_Spatial()

    # Combinaciones posibles de años y trimestres para la interpolación
    combinaciones <- purrr::transpose(
        purrr::cross2(.x = unique(simulation_dates$year),
                      .y = unique(simulation_dates$season)))

    # Interpolación de acumulados de lluvia
    datos_interpolados <- purrr::pmap_dfr(
        .l = combinaciones,
        .f = function(year, season) {

            interpolacion.valores.iniciales.sp <- datos.estacionales %>%
                dplyr::filter(year == !!year & season == !!season) %>%
                sf::as_Spatial()

            # Valores interpolados de ocurrencia
            datos_interpolados_prcp <-
                automap::autoKrige(seasonal_prcp~longitude+latitude,
                                   interpolacion.valores.iniciales.sp,
                                   simulation_points.sp,
                                   debug.level = 0) %>%
                sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
                dplyr::mutate(var1.pred = if_else(var1.pred < 0, 0, var1.pred)) %>%
                dplyr::select(seasonal_prcp = var1.pred)


            # Valores interpolados de temperatura máxima
            datos_interpolados_tmax <-
                automap::autoKrige(seasonal_tmax~longitude+latitude,
                                   interpolacion.valores.iniciales.sp,
                                   simulation_points.sp,
                                   debug.level = 0) %>%
                sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
                dplyr::select(seasonal_tmax = var1.pred)


            # Valores interpolados de temperatura mínima
            datos_interpolados_tmin <-
                automap::autoKrige(seasonal_tmin~longitude+latitude,
                                   interpolacion.valores.iniciales.sp,
                                   simulation_points.sp,
                                   debug.level = 0) %>%
                sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
                dplyr::select(seasonal_tmin = var1.pred)


            # Crear data frame con los valores interpolados
            # Este df hay que unirlo a la grilla de simulación para cada día de la simulacion
            # para la union tenemos las coordenadas, anos y estaciones.
            datos_interpolados_prcp_tmax_tmin <- simulation_points %>%
                dplyr::select(longitude, latitude) %>%
                sf::st_join(datos_interpolados_prcp) %>%
                sf::st_join(datos_interpolados_tmax) %>%
                sf::st_join(datos_interpolados_tmin) %>%
                dplyr::mutate(year = !!year, season = !!season) %>%
                dplyr::select(year, season, seasonal_prcp, seasonal_tmax, seasonal_tmin, longitude, latitude) %>%
                sf::st_drop_geometry() %>% tibble::as_tibble()

            return(datos_interpolados_prcp_tmax_tmin)
        }
    )

    return (datos_interpolados %>%
                sf::st_as_sf(coords = c("longitude", "latitude")) %>%
                sf::st_set_crs(sf::st_crs(simulation_points)))

}


# Generacion de campos gaussianos para la alteración de los valores iniciales,
# haciendo uso de los datos utilizados para realizar el ajuste
generate_variograms_for_initial_values <- function(model, simulation_points, seed, month, day) {

    # Creación de la matriz de simulación
    simulation_matrix <- model$models_data %>%
        dplyr::mutate(year = lubridate::year(date),
                      month = lubridate::month(date),
                      day = lubridate::day(date),
                      time = as.numeric(date)/1000) %>%
        dplyr::left_join(model$stations, by = "station_id") %>%
        dplyr::select(station_id, date, year, month, day,
                      season, prcp_occ, tmax, tmin, type_day,
                      longitude, latitude, geometry) %>%
        sf::st_as_sf(crs = sf::st_crs(model$stations))

    # Creación de la matriz de distancias
    distance_matrix <- gamwgen:::make_distance_matrix(model$stations)

    # Generación de variogramas
    variograms_for_initial_values <-
        gamwgen:::setting_variograms_for_initial_values(
            simulation_matrix = simulation_matrix,
            distance_matrix = distance_matrix,
            grid = simulation_points %>%
                dplyr::select(longitude, latitude) %>%
                sf::st_drop_geometry() %>%
                tibble::as_tibble(),
            seed = seed,
            init_values_month = month,
            init_values_day = day)

    return (variograms_for_initial_values)

}

interpolate_start_climatology <- function(model, simulation_points, seed, month, day) {

    # OBS:
    # Al interpolar los datos para el día previo al día de inicio de la simulación, no se
    # tiene en cuenta si el usario optó por usar ruidos correlacionados espacialmente o no.
    # Esto es correcto, es el funcionamiento previsto y no un error!

    # Generación de valores inciales para el variograma
    variograms_for_initial_values <-
        gamwgen:::generate_variograms_for_initial_values(model, simulation_points, seed, month, day)

    # Matriz de valores a interpolar
    data_to_be_interpolated <- model$start_climatology %>%
        dplyr::filter(month == !!month, day == !!day) %>%
        dplyr::left_join(model$stations, by = "station_id") %>%
        sf::st_as_sf(crs = sf::st_crs(model$stations))

    # Valores interpolados de ocurrencia
    prcp_occ_interpolation <-
        automap::autoKrige(prcp_occ~1,
                           data_to_be_interpolated %>% sf::as_Spatial(),
                           simulation_points %>% sf::as_Spatial(),
                           debug.level = 0) %>%
        sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
        dplyr::mutate(noise = variograms_for_initial_values$random_fields$prcp,
                      prcp_occ = as.numeric(var1.pred + noise > 0))

    # Valores interpolados de temperatura maxima
    tmax_interpolation <-
        automap::autoKrige(tmax~1,
                           data_to_be_interpolated %>% sf::as_Spatial(),
                           simulation_points %>% sf::as_Spatial(),
                           debug.level = 0) %>%
        sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
        dplyr::mutate(noise = variograms_for_initial_values$random_fields$tmax,
                      tmax = var1.pred + noise)

    # Valores interpolados de temperatura mínima
    tmin_interpolation <-
        automap::autoKrige(tmin~1,
                           data_to_be_interpolated %>% sf::as_Spatial(),
                           simulation_points %>% sf::as_Spatial(),
                           debug.level = 0) %>%
        sf::st_as_sf(x = .$krige_output, crs = sf::st_crs(simulation_points)) %>%
        dplyr::mutate(noise = variograms_for_initial_values$random_fields$tmin,
                      tmin = var1.pred + noise)

    return (prcp_occ_interpolation %>% dplyr::select(prcp_occ) %>%
                sf::st_join(tmax_interpolation %>% dplyr::select(tmax)) %>%
                sf::st_join(tmin_interpolation %>% dplyr::select(tmin)))
}


# krige_covariate <- function(model, station_coordinates_grid, simulation_coordinates_grid, covariate) {
#     variogram <- fields::vgram(loc = station_coordinates_grid, y = covariate,
#                                breaks = seq(0, max(model$distance_matrix),
#                                             by = max(model$distance_matrix) / 30))
#     variogram <- data.frame(y = variogram$stats[2, ], x = variogram$centers)
#
#     krigged_coefs <- tryCatch({
#         # variogram_fit <- with(variogram,
#         #                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#         #                                  start = list(rho = max(y, na.rm = TRUE), theta = max(dist_values)/2),
#         #                                  control = list(maxiter = 500, minFactor = 1/10000, warnOnly = F)))
#         variogram_fit <- with(variogram,
#                               stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                                          algorithm = "port",
#                                          start = list(rho = max(y, na.rm = TRUE), theta = max(model$distance_matrix)/2),
#                                          lower = list(rho = 0, theta = 0),
#                                          upper = list(rho = max(y, na.rm = TRUE) * 1.5, theta = max(model$distance_matrix) * 1.5)))
#         krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F,
#                                    sigma2 = 0, rho = coefficients(variogram_fit)['rho'], theta = coefficients(variogram_fit)['theta'])
#         fields::predict.Krig(krig_model, simulation_coordinates_grid)
#
#     }, error = function(...) {
#         fields::predict.Krig(fields::Krig(station_coordinates_grid, covariate, give.warnings = F), simulation_coordinates_grid)
#     })
# }


krige_covariate_automap <- function(model, stations_locations, simulation_locations, covariate) {
    covariate <- SpatialPointsDataFrame(sp::coordinates(stations_locations),
                                        data.frame(value = covariate),
                                        proj4string = stations_locations@proj4string)
    automap::autoKrige(value ~ 1, input_data = covariate, new_data = simulation_locations, debug.level = 0)$krige_output$var1.pred
}

idw_covariate <- function(model, stations_locations, simulation_locations, covariate) {
    drop(model$idw_weights %*% covariate)
}


## TODO: cleanup.

# estimate model coefficients on grid using ordinary kriging (OK) occurrence

# dist_function <- ifelse(is_projected, 'rdist', 'rdist.earth')
#
#
#         krig_covariate <- function(covariate) {
#             # krig_model <- fields::Krig(coordinates(model$stations), covariate, give.warnings = F, Distance = dist_function)
#             krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F)
#             fields::predict.Krig(krig_model, simulation_coordinates_grid)
#         }
#
#         krig_gstat <- function(covariate) {
#             covariate <- data.frame(covariate, coordinates(model$stations))
#             colnames(covariate) <- c('value', 'x', 'y')
#             coordinates(covariate) <- ~ x + y
#             covariate@proj4string <- model$stations@proj4string
#
#             max_range <- max(model$distance_matrix)
#
#             vgm_model <- c('Exp', 'Gau', 'Log')
#
#             result <- NULL
#             model_index <- 0
#
#             krige_internal <- function(vgm_model) {
#                 vario <- gstat::variogram(value ~ 1, data = covariate, cutoff = max_range, width = max_range / 20)
#                 # vario <- gstat::variogram(value ~ 1, data = covariate)
#                 vario_fit <- gstat::fit.variogram(vario, gstat::vgm(model = vgm_model[model_index],
#                                                                     psill = mean(vario$gamma, na.rm = T),
#                                                                     range = max_range / 10))
#                 krige_fit <- gstat::krige(value ~ 1, covariate, krige_locations, model = vario_fit, debug.level = 0)
#                 krige_fit
#             }
#
#             while(is.null(result) && model_index < length(vgm_model)) {
#                 model_index <- model_index + 1
#                 result <- tryCatch(krige_internal(vgm_model), warning = function(w) {
#                     NULL
#                 })
#             }
#
#             if(is.null(result)) {
#                 warning('Failed to build variogram')
#                 return(NA)
#             } else {
#                 covariate_idx <- apply(spDists(covariate, result), 1, which.min)
#                 return(result[covariate_idx, ]$var1.pred)
#             }
#         }
#
#         # krig_automap <- function(covariate) {
#         #     covariate <- data.frame(covariate, coordinates(model$stations))
#         #     colnames(covariate) <- c('value', 'x', 'y')
#         #     coordinates(covariate) <- ~ x + y
#         #     covariate@proj4string <- model$stations@proj4string
#         #     automap::autoKrige(value ~ 1, covariate, krige_locations)
#         # }
#
#         spline_covariate <- function(covariate) {
#             spline_model <- fields::Tps(coordinates(model$stations), covariate, give.warnings = F, lon.lat = !is_projected, miles = F)
#             fields::predict.Krig(spline_model, simulation_coordinates)
#         }
#
#         dist_values <- do.call(c, sapply(1:nrow(model$distance_matrix), function(r) model$distance_matrix[r, (r:nrow(model$distance_matrix))]))
#         dist_values <- dist_values[dist_values > 0]
#         #
#         # covariate_sp <- station_sp_gk
#         # covariate_sp@data$value <- covariate
#         # automap::autofitVariogram(value ~ 1, covariate_sp)
#
#         krig_covariate_2 <- function(covariate) {
#             # variogram <- fields::vgram(loc = station_coordinates_grid, d = dist_values, y = covariate,
#             #                            breaks = seq(0, max(dist_values), by = max(dist_values) / 20))
#             variogram <- fields::vgram(loc = station_coordinates_grid, y = covariate,
#                                        breaks = seq(0, max(model$distance_matrix), by = max(model$distance_matrix) / 30))
#             variogram <- data.frame(y = variogram$stats[2, ], x = variogram$centers)
#
#             krigged_coefs <- tryCatch({
#                 # variogram_fit <- with(variogram,
#                 #                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                 #                                  start = list(rho = max(y, na.rm = TRUE), theta = max(dist_values)/2),
#                 #                                  control = list(maxiter = 500, minFactor = 1/10000, warnOnly = F)))
#                 variogram_fit <- with(variogram,
#                                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                                                  algorithm = "port",
#                                                  start = list(rho = max(y, na.rm = TRUE), theta = max(model$distance_matrix)/2),
#                                                  lower = list(rho = 0, theta = 0),
#                                                  upper = list(rho = max(y, na.rm = TRUE) * 1.5, theta = max(model$distance_matrix) * 1.5)))
#                 krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F,
#                                            sigma2 = 0, rho = coefficients(variogram_fit)['rho'], theta = coefficients(variogram_fit)['theta'])
#                 fields::predict.Krig(krig_model, simulation_coordinates_grid)
#
#             }, error = function(...) {
#                 fields::predict.Krig(fields::Krig(station_coordinates_grid, covariate, give.warnings = F), simulation_coordinates_grid)
#             })
#
#         }
