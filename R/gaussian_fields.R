
# Definicion de funcion para la generacion de campos gaussianos para alterar los valores iniciales
setting_variograms_for_initial_values <- function(simulation_matrix, distance_matrix, grid, seed,
                                                  init_values_month = 12, init_values_day = 31) {

    ## OBS:
    # Estamos trabajando con coordenadas métricas, sin embargo,
    # el paquete RandomFields trabaja mejor en metros, por la tanto,
    # pasamos distance_matrix y grid a kilometros!!

    # aquí pasamos distance_matrix a kilómetros!!
    distance_matrix <- distance_matrix / 1000
    # y aquí pasamos grid a kilómetros!!
    grid <- grid %>%
        dplyr::mutate(longitude = longitude / 1000,
                      latitude  = latitude / 1000)

    # simulation matrix debe ser un tibble (no un sf)
    simulation_matrix <- simulation_matrix %>%
        tibble::as_tibble()

    # Matriz de valores iniciales para inicializar la simulacion

    prcp_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(month == init_values_month, day == init_values_day) %>%
        dplyr::select(station_id, date, prcp_occ) %>%
        tidyr::spread(key = station_id, value = prcp_occ) %>%
        dplyr::select(colnames(distance_matrix)) %>%
        dplyr::mutate_if(is.numeric, as.factor)

    prcp_intial_values_interpolation <- missMDA::imputeMCA(
        don = prcp_intial_values_interpolation, ncp = 5)$completeObs %>%
        dplyr::mutate_if( is.factor, as.character) %>%
        dplyr::mutate_if( is.character, as.numeric)

    tmax_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(month == init_values_month, day == init_values_day) %>%
        dplyr::select(station_id, date, tmax) %>%
        tidyr::spread(key = station_id, value = tmax) %>%
        dplyr::select(colnames(distance_matrix))

    tmax_intial_values_interpolation <- missMDA::imputePCA(
        X = tmax_intial_values_interpolation)$completeObs

    tmin_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(month == init_values_month, day == init_values_day) %>%
        dplyr::select(station_id, date, tmin) %>%
        tidyr::spread(key = station_id, value = tmin) %>%
        dplyr::select(colnames(distance_matrix))

    tmin_intial_values_interpolation <- missMDA::imputePCA(
        X = tmin_intial_values_interpolation)$completeObs


    prcp_cor   <- stats::cor(prcp_intial_values_interpolation, use = "pairwise.complete")
    prcp_vario <- stats::var(prcp_intial_values_interpolation, use = "pairwise.complete") * (1 - prcp_cor)
    tmax_vario <- stats::cov(tmax_intial_values_interpolation, use = "pairwise.complete.obs")
    tmin_vario <- stats::cov(tmin_intial_values_interpolation, use = "pairwise.complete.obs")

    prcp_params <- stats::optim(par = c(0.01, 1, max(distance_matrix)),
                                f = gamwgen:::partially_apply_LS(prcp_vario, distance_matrix))$par
    prcp_params[prcp_params < 0] <- 0
    prcp_params <- c(0, prcp_params[2], prcp_params[3])

    sill_initial_value_tmax <- mean(var(tmax_intial_values_interpolation, na.rm = T))
    tmax_params <- stats::optim(par = c(sill_initial_value_tmax, max(distance_matrix)),
                                fn = gamwgen:::partially_apply_LS(tmax_vario, distance_matrix, base_p = c(0)))$par
    tmax_params <- c(0, tmax_params)

    sill_initial_value_tmin <- mean(var(tmin_intial_values_interpolation, na.rm = T))
    tmin_params <- stats::optim(par = c(sill_initial_value_tmin, max(distance_matrix)),
                                fn = gamwgen:::partially_apply_LS(tmin_vario, distance_matrix, base_p = c(0)))$par
    tmin_params <- c(0, tmin_params)

    variogram_parameters <- list(prcp = prcp_params,
                                 tmax = tmax_params,
                                 tmin = tmin_params)

    models_parameters <- list(prcp = RandomFields::RMexp(var = variogram_parameters$prcp[2],
                                                         scale = variogram_parameters$prcp[3]),
                              tmax = RandomFields::RMexp(var = variogram_parameters$tmax[2],
                                                         scale = variogram_parameters$tmax[3]),
                              tmin = RandomFields::RMexp(var = variogram_parameters$tmin[2],
                                                         scale = variogram_parameters$tmin[3]))

    random_fields <- list(prcp = RandomFields::RFsimulate(models_parameters$prcp, x = grid, grid = F,
                                                          seed = seed, coord_sys="cartesian", coordunits="km"),
                          tmax = RandomFields::RFsimulate(models_parameters$tmax, x = grid, grid = F,
                                                          seed = seed, coord_sys="cartesian", coordunits="km"),
                          tmin = RandomFields::RFsimulate(models_parameters$tmin, x = grid, grid = F,
                                                          seed = seed, coord_sys="cartesian", coordunits="km"))

    return(list(variogram = variogram_parameters,
                models = models_parameters,
                random_fields = random_fields))

}
